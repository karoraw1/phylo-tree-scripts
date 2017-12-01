# -*- coding: utf-8 -*-
# define dictionary of species and accession relationships
from collections import OrderedDict
import urllib2, time, sys, os
import subprocess as sp
from Bio import SeqIO

skips = [False, False, False, False]
if "--skipDL" in sys.argv:
	skips[0] = True
if "--skipAL1" in sys.argv:
	skips[1] = True
if "--skipMF" in sys.argv:
	skips[2] = True
if "--skipAL2" in sys.argv:
	skips[3] = True


informative_subset = [37, 127, 223, 68, 54, 334, 97, 65, 461, 373, 316, 415, 431, 3, 1, 352, 470, 120, 271, 37, 189,
                      209, 238, 191, 268, 123, 88, 463]

type_lists = OrderedDict()
type_lists["Escherichia coli"] = ["J01859"]
type_lists["Methanosarcina mazei"] = ["AB065295"]
type_lists["Desulfobulbus japonicus"] = ["AB110550"]
type_lists["Desulfobulbus rhabdoformis"] = ["AB546247"]
type_lists["Desulfotalea psychrophila"] = ["AF099062"]
type_lists["Desulfobulbus mediterraneus"] = ["AF354663"]
type_lists["Desulfobacterium corrodens"] = ["AY274450"]
type_lists["Desulfobulbus propionicus"] = ["AY548789"]
type_lists["Desulfovibrio marinus"] = ["DQ365924"]
type_lists["Desulfurivibrio alkaliphilus"] = ["EF422413"]
type_lists["Desulfobulbus alkaliphilus"] = ["HM750216"]
type_lists["Desulfobulbus elongatus"] = ["X95180"]
type_lists["Desulfocapsa sulfoexigens"] = ["Y13672"]
type_lists["Electrothrix marina"] = ["JX091065", "JX091056", "JX091026", "KP265606", "KR912340", 
                                     "KR912341", "KR912342"]
type_lists["Electrothrix communis"] = ["KJ021898", "KJ021897", "KJ021896", "KJ021895", "KJ021894", 
									   "JX091070", "JX091067", "JX091062", "HG004417", "HG004416",
									   "JX091054", "JX091052", "JX091041", "JX091028", "JX091025", 
									   "HG004415", "HG004414", "HG004413", "HG004412", "HG004411", 
									   "HG004410", "HG004409", "HG004408", "HG004407", "HG004405", 
									   "HG004425", "HG004420", "HG004419", "KR912339", "KR912343", 
									   "KR912344", "KR912345", "KR912346", "KR912347", "KR912348"]
type_lists["Electronema palustris"] = ["FQ658891", "FQ658831", "GU208270"]
type_lists["Electrothrix japonica"] = ["KJ562733", "KP265514", "JX091073", "JX091064", "JX091057", 
                                       "HG004406", "HG004418", "GQ249497", "JF268345", "JF268368", 
                                       "JF268348", "KR912349"]
type_lists["Electrothrix aarhusiensis"] = ["HG004404", "KR912338"]


# download files

main_file = "type_strains.fasta"
url_template = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=%s&rettype=fasta&retmode=text"

sp_lookup = {}

if not skips[0]:
	wipe = open(main_file, "w")
	wipe.close()

for species, id_list in type_lists.items():
	for id in id_list:
		sp_lookup[id] = species
		if not skips[0]:
			open(main_file, "a").write(urllib2.urlopen(url_template % id).read())
			time.sleep(1.0/3)
			sys.stdout.write("Fetched %s...\n" % id)
    		sys.stdout.flush()

def perform_alignment(in_name, report_name):
	in_prefix = in_name.split(".")[0]
	aln_cmd_1 = "cmalign --outformat Pfam -o " + in_prefix + ".align.sto "
	hmm_model = "bacteria_ssu.cm "
	aln_cmd = aln_cmd_1 + hmm_model + in_name + " > " + report_name
	align_16S = sp.Popen(aln_cmd, cwd = os.getcwd(), shell = True, stderr=sp.PIPE, stdout=sp.PIPE)
	stdout_stderr = align_16S.communicate()
	return stdout_stderr

if not skips[1]:
	_ = perform_alignment(main_file, "alignment_report.txt")


with open("alignment_report.txt", "r") as arh:
    aln_result = arh.read()
aln_lines = aln_result.split("\n")[16:-3]

not_included = type_lists.keys()

good_types = OrderedDict()


for a_line in aln_lines:
	a_tabs = a_line.split()
	accs = a_tabs[1].split(".")[0]
	this_sp = sp_lookup[accs]
	a_start = int(a_tabs[3])
	an_end = int(a_tabs[4])
	total_len = int(a_tabs[2])

	if (a_start > 793) or (an_end < 524):
		the_overlap = 0.0
	elif (a_start < 524) and (an_end > 793):
		the_overlap = 1.0
	else:
		front_gap, end_gap, front_trim, end_trim = 0, 0, 0, 0
		if a_start > 524:
			front_gap = a_start - 524
		if an_end < 793:
			end_gap = 793 - an_end
		the_overlap = round( ((270. - front_gap - end_gap) / 270.), 2)
	if the_overlap > 0.95:
		if a_start < 524:
			front_trim = 524 - a_start
		if an_end > 793:
			end_trim = an_end - 793
		good_types[accs] = (front_trim, end_trim)
	if this_sp in not_included:
		not_included.remove(this_sp)
    
print not_included

good_accs = good_types.keys()
sys.stdout.write("%r sequences retained of %r\n" % (len(good_accs), len(sp_lookup.keys())))
sys.stdout.flush()

# Create a file containing type strains and new tags
mixed_seqs_f = "type_and_new_str.fasta"

if not skips[2]:
	ches_seqs_f = "chesapeake_desulf.fasta"
	new_sequences = SeqIO.parse(open(ches_seqs_f, "r"), 'fasta')
	type_sequences = SeqIO.parse(open(main_file, "r"), 'fasta')

	# write all the chesapeake tags into a new mixed file
	with open(mixed_seqs_f, "w") as of_:
	    for new_seq in new_sequences:
	    	name, sequence = new_seq.id, str(new_seq.seq)
	    	if int(name) in informative_subset:
	    		# modify the names to add "Seq"
	    		name = ">Seq_" + name
	    		of_.write(name+"\n")
	    		of_.write(sequence+"\n")

	# go through type strain file     
	# identify the members of the good hit club
	    for type_seq in type_sequences:
	    	t_name, t_seq = type_seq.id, str(type_seq.seq)
	    	match_name = t_name.split(".")[0]
	    	if match_name in good_accs:
	    		# trim from front and back 
	    		clean_seq = t_seq[good_types[match_name][0]: -good_types[match_name][1]]
		    	# modify the header so it just has the accesion and the name
		    	clean_name = ">"+match_name+" "+sp_lookup[match_name]
		    	# copy into the mixed file 
		    	of_.write(clean_name+"\n")
		    	of_.write(clean_seq+"\n")

# redo the alignment
if not skips[3]:
    _ = perform_alignment(mixed_seqs_f, "trimmed_aln_rep.txt")

# parse the new alignment report to verify trimming
	#TODO

# convert the new alignment into a fasta

aln_root = "type_and_new_str.align"
converted_file = aln_root + '.fasta'
if not os.path.exists(converted_file):
	seqmag_cmd = 'seqmagick convert ' + aln_root+ '.sto '+ aln_root + '.fasta'
	convert = sp.Popen(seqmag_cmd, cwd = os.getcwd(), shell = True, stderr=sp.PIPE, stdout=sp.PIPE)
	convert.communicate()
else:
	sys.stdout.write("Converted file detected!")
	sys.stdout.flush()

# create ML distance data
tree_file = "RAxML_parsimonyTree.dist"
if not os.path.exists(tree_file):
	raxml_cmd = "./raxml -T 4 -f x -p 12345 -s " + converted_file + ' -m GTRGAMMA -n dist'
	dist = sp.Popen(raxml_cmd, shell = True, stderr=sp.PIPE, stdout=sp.PIPE)
	dist.communicate()

new_tree_file = "newick_tree_modified.dist"

with open(tree_file, "r") as old_tree:
	full_tree = old_tree.read()

for each_accs in good_accs:
	new_label = (each_accs + "_" + sp_lookup[each_accs]).replace(" ", "-")
	full_tree = full_tree.replace(each_accs, new_label)

with open(new_tree_file, "w") as new_tree_h:
	new_tree_h.write(full_tree)

