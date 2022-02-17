#!/home/rs2474/anaconda3/envs/general/bin/python

### Use python to output SLiM file for running Buffering Base models for fitting empirical data by varying RC and RS

import sys
import subprocess
import time
import argparse
import pyfaidx as px
from collections import Counter


def wrt2buff(buff_file, conserv_aa_ratio, sel_aa_ratio, conserv_coeff, pop_sample, sel_coeff, phase_len=12500, neutral_burnIn_len=20000, wolbachia_non_burnIn_len=25000, wolbachia_first=True):
    """
        wrt2buff: 
            write the buffering-selection regime slim model
        
        Arguments:
            buff_file (str): the slim file to output
            conserv_aa_ratio (float): the proportion of conserved AAs
            sel_aa_ratio (float): the proportion of selected AAs 
            conserv_coeff (float): the selection coefficient (<0) for non-synonymous mutation at conserved sites
            pop_sample (int): the number of diploid individuals for population sample
            sel_coeff (float): the selection coefficient (absolute value) for alternating wolbachia infection phases
            phase_len (int): the length of each wolbachia/non-wolbachia infection period
            neutral_burnIn_len (int): the length of neutral burn-in period for accumulating neutral polymorphisms
            wolbachia_non_burnIn_len (int): the total length of alternating wolbachia infection phases
            wolbachia_first (bool): wolbachia infection phase kicks in first (True) or non-wolbachia infection phase kicks in first (False)
            
        Return:
            None
    """
    num_phases = int(wolbachia_non_burnIn_len) // phase_len
    ## register the 1st, 3rd, 5th, etc. phases in Wolbachia alternating stages
    first_phases = ','.join(list(map(lambda x: '{}:{}'.format(neutral_burnIn_len+phase_len*x+1,neutral_burnIn_len+x*phase_len+phase_len), range(0,num_phases-1,2))))
    ## register the 2nd, 4th, 6th, etc. phases in Wolbachia alternating stages
    second_phases= ','.join(list(map(lambda x: '{}:{}'.format(neutral_burnIn_len+phase_len*x+1, neutral_burnIn_len+x*phase_len+phase_len), range(1,num_phases,2))))
    with open(buff_file, 'w') as f:
        f.write("// function to read in the miyata matrix\n")
        f.write("function (f)readFloatTable(s$ path){\n")
        f.write("\tif (!fileExists(path))\n")
        f.write("\t\tstop(\"readFloatTable(): File not found at path \" + path);\n")
        f.write("\tl = readFile(path);\n")
        f.write('\tm = sapply(l, "asFloat(strsplit(applyValue, sep=\''+r'\t'+'\'));", simplify="matrix");'+'\n')
        f.write("\treturn t(m);\n")
        f.write("}\n\n\n")
        
        
        f.write("// check if there's an early stop codon or no stop codon in the two genomes of an individual\n")
        f.write("function (f)checkStopCodon(object<Individual> ind) {\n")
        f.write("\t// check both genomes of the individual to see if there's any stop codon in any of them\n")
        f.write("\tfor (g in ind.genomes) {\n")
        f.write("\t\t// get the sequence before stop codon\n")
        f.write("\t\tnucSeq = g.nucleotides(0, SeqLen-4);\n")
        f.write("\t\tseq2Codons = nucleotidesToCodons(nucSeq);\n")
        f.write("\t\t// get the stop codon seq\n")
        f.write("\t\tstopSeq = nucleotidesToCodons(g.nucleotides(SeqLen-3, SeqLen-1));\n")
        f.write("\t\t// early stop codons OR no stop codon OR no start codon -> lethal\n")
        f.write("\t\tif (any(match(seq2Codons, STOP)!=-1) | !(match(stopSeq, STOP)>=0) | substr(nucSeq,0,2)!='ATG')\n")
        f.write("\t\t\treturn 0.0;\n")
        f.write("\t\t}\n")
        f.write("\t\t// if no early stop codon was detected and this function got returned,\n")
        f.write("\t\t// individual successfully survive for potential further fitness modification\n")
        f.write("\t\treturn 1.0;\n")
        f.write("}\n\n\n")
        
        
        f.write("// set up a bam simulation model\n")
        f.write("initialize() {\n")
        f.write("\t// define some constants to be used later\n")
        f.write("\tdefineConstant(\"TAA\", nucleotidesToCodons(\"TAA\"));\n")
        f.write("\tdefineConstant(\"TAG\", nucleotidesToCodons(\"TAG\"));\n")
        f.write("\tdefineConstant(\"TGA\", nucleotidesToCodons(\"TGA\"));\n")
        f.write("\tdefineConstant(\"STOP\", c(TAA, TAG, TGA));\n")
        f.write("\tdefineConstant(\"NONSTOP\", (0:63)[match(0:63, STOP) < 0]);\n\n")
        f.write("\t// selection coefficients predefined\n")
        f.write("\tdefineConstant(\"DelSelCoeff\", -{});\n".format(sel_coeff))
        f.write("\tdefineConstant(\"BenSelCoeff\", {});\n".format(sel_coeff))
        f.write("\tdefineConstant(\"ConservCoeff\", {});\n\n".format(conserv_coeff))
        
        f.write("\t// initialize a nucleotide-based SLiM model\n")
        f.write("\tinitializeSLiMOptions(nucleotideBased=T);\n")
        f.write("\tinitializeSex(\"A\");\n")
        f.write("\t// read in the reference nucleotide sequence\n")
        f.write("\tL = initializeAncestralNucleotides('melsim_ancestral.fasta');\n")
        f.write("\tdefineConstant(\"SeqLen\", L);\n")
        f.write("\t// get the 85% of the AAs that are identical between species\n")
        f.write("\tsame_sites_line = readFile('mel-sim_conserv_sites.txt');\n")
        f.write('\tsameAAs = sapply(same_sites_line, "asInteger(strsplit(applyValue, sep=\''+r'\t'+'\'));");\n')
        f.write("\tcat(\"SameSites: \");\n\tcatn(sameAAs);\n")
        f.write('\tdiffAAs = setDifference(0:integerDiv(L,3), sameAAs); // exclude start and stop codons already\n')
        
        f.write("\t// randomly select the AAs to be presumably conserved, neutral and under selection\n")
        f.write('\tnum_conserv_AAs = asInteger(integerDiv(L,3)*{});\n'.format(conserv_aa_ratio))
        f.write('\tnum_sel_AAs = asInteger(integerDiv(L,3)*{});\n'.format(sel_aa_ratio))
        
        f.write('\tconservAAs = sample(sameAAs, num_conserv_AAs);\n')
        f.write("\t// convert the conserved AAs' positions to codons' nucleotide positions\n")
        f.write("\tdefineConstant(\"ConservSites\", c(conservAAs*3, conservAAs*3+1, conservAAs*3+2));\n")
        f.write("\t// print out the conservSites\n")
        f.write("\tcat(\"ConservSites: \");\n\tcatn(ConservSites);\n")
        
        f.write("\t// sample the AAs under selection\n")
        f.write('\tsampled_selAAs = sample(setDifference(0:integerDiv(L,3), conservAAs), num_sel_AAs); // exclude start and stop codons\n')
        f.write("cat(\"SelectedSites: \");\n\tcatn(c(sampled_selAAs*3, sampled_selAAs*3+1, sampled_selAAs*3+2));\n")
        
        f.write("\t// set the AAs to be neutral\n")
        f.write('\tneutralAAs = setDifference(c(0:integerDiv(L,3)), c(conservAAs, sampled_selAAs));\n')
        f.write("\tdefineConstant(\"NeutralSites\", c(neutralAAs*3, neutralAAs*3+1, neutralAAs*3+2));\n")
        f.write("\tcat(\"NeutralSites: \");\n\tcatn(NeutralSites);\n")
        f.write("\t// read in the 64x64 codon-based Miyata distance matrix\n")
        f.write("\tdefineConstant(\"M\", readFloatTable(\"codon_miyata_mat.txt\"));\n")
        f.write("\t// initialize the neutral mutations\n")
        f.write("\tinitializeMutationTypeNuc(\"m1\", 0.5, \"f\", 0.0);\n")
        f.write("\t// initialize the mutation rate with mmKimura to characterize transition vs transversion rate\n")
        f.write("\tinitializeGenomicElementType(\"g1\", m1, 1.0, mmKimura(1.88e-6, 0.47e-6));\n")
        f.write("\t// initialize the genomic segment and uniform recombination rate\n")
        f.write("\tinitializeGenomicElement(g1, 0, L-1);\n")
        f.write("\tinitializeRecombinationRate(1e-5);\n")
        f.write("}\n\n\n")
        
        
        f.write("// define a fitness(NULL) callback that makes early stop codons lethal\n")
        f.write("fitness(NULL) {\n")
        f.write("\treturn checkStopCodon(individual);\n")
        f.write("}\n\n\n")
        
        
        f.write("// define a mutation callback to control mutation generation and set the non-synonymous mutations at conserved sites as deleterious mutations during the neutral burn-in period\n")
        f.write("1:%d mutation(m1) {\n" % neutral_burnIn_len)
        f.write("\tsim.setValue(\"TotalMuts\", sim.getValue(\"TotalMuts\")+1);\n")
        f.write("\tsim_mut_aa_pos = integerDiv(sim.mutations.position, 3);\n")
        f.write("\tsim_mut_nuc_pos = c(sim_mut_aa_pos*3, sim_mut_aa_pos*3+1, sim_mut_aa_pos*3+2);\n")
        f.write("\t// if the mutation occurs at the codon already bearing a mutation and the site isn't conserved/neutral, block it\n")
        f.write("\tif (match(mut.position, sim_mut_nuc_pos)!=-1 & match(mut.position, ConservSites)==-1 & match(mut.position, NeutralSites)==-1){\n")
        f.write("\t\tsim.setValue(\"RejMuts\", sim.getValue(\"RejMuts\")+1);\n")
        f.write("\t\treturn F;\n")
        f.write("\t}\n")
        
        f.write("\t// get the original parental genome before mutation\n")
        f.write("\toriginal_genome = genome.nucleotides();\n")
        f.write("\t// get the new genome with mutation\n")
        f.write("\tmut_genome = substr(original_genome, 0, mut.position-1)+mut.nucleotide+substr(original_genome, mut.position+1);\n")
        f.write("\t// get the mut codon and original codon\n")
        f.write("\tmut_codon_pos = integerDiv(mut.position, 3);\n")
        f.write("\toriginal_codon_nuc = substr(original_genome, mut_codon_pos*3, mut_codon_pos*3+2);\n")
        f.write("\tmut_codon_nuc = substr(mut_genome, mut_codon_pos*3, mut_codon_pos*3+2);\n")
        f.write("\t// get the mutated and original codons and AAs respectively\n")
        f.write("\toriginal_codon = nucleotidesToCodons(original_codon_nuc);\n")
        f.write("\tmut_codon = nucleotidesToCodons(mut_codon_nuc);\n")
        f.write("\toriginal_AA = codonsToAminoAcids(original_codon);\n")
        f.write("\tmut_AA = codonsToAminoAcids(mut_codon);\n")
        f.write("\t// check if non-syn mutation and the site is conserved \n")
        f.write("\tif ((original_AA!=mut_AA) & match(mut.position, ConservSites)!=-1) {\n")
        f.write("\t\tmut.setSelectionCoeff(ConservCoeff);\n")
        f.write("\t}\n")
        f.write("\treturn T;\n")
        f.write("}\n\n\n")
        
        
        f.write("// define a mutation callback during Wolbachia infection phase to control mutation generation and reset the new mutation's selection coefficient based on the Miyata distance between the mutated amino acid and the original\n")
        ## if wolbachia-infection phase comes first, then this section is defined to be the s1 section
        if (wolbachia_first):
            f.write("s1 %d:%d mutation(m1) {\n" % (neutral_burnIn_len+1, neutral_burnIn_len+phase_len))
        ## if wolbachia-lost phase comes first, then this section is defined to be the s2 section
        else:
            f.write("s2 %d:%d mutation(m1) {\n" % (neutral_burnIn_len+phase_len+1, neutral_burnIn_len+2*phase_len))
        f.write("\tsim.setValue(\"TotalMuts\", sim.getValue(\"TotalMuts\")+1);\n")
        f.write("\tsim_mut_aa_pos = integerDiv(sim.mutations.position, 3);\n")
        f.write("\tsim_mut_nuc_pos = c(sim_mut_aa_pos*3, sim_mut_aa_pos*3+1, sim_mut_aa_pos*3+2);\n")
        f.write("\t// if the mutation occurs at the codon already bearing a mutation and the site isn't conserved/neutral, block it\n")
        f.write("\tif (match(mut.position, sim_mut_nuc_pos)!=-1 & match(mut.position, ConservSites)==-1 & match(mut.position, NeutralSites)==-1){\n")
        f.write("\t\tsim.setValue(\"RejMuts\", sim.getValue(\"RejMuts\")+1);\n")
        f.write("\t\treturn F;\n")
        f.write("\t}\n")
        
        f.write("\t// if the mutation occurs at the neutral codon, return it with s = 0\n")
        f.write("\tif (match(mut.position, NeutralSites)!=-1){\n")
        f.write("\t\tmut.setSelectionCoeff(0.0);\n")
        f.write("\t\treturn T;\n")
        f.write("\t}\n")
        
        f.write("\t// get the original parental genome before mutation\n")
        f.write("\toriginal_genome = genome.nucleotides();\n")
        f.write("\t// get the new genome with mutation\n")
        f.write("\tmut_genome = substr(original_genome, 0, mut.position-1)+mut.nucleotide+substr(original_genome, mut.position+1);\n")
        f.write("\t// get the mutated and original codons\n")
        f.write("\toriginal_codons = nucleotidesToCodons(original_genome);\n")
        f.write("\tmut_codons = nucleotidesToCodons(mut_genome);\n")
        f.write("\toriginal_AAs = codonsToAminoAcids(original_codons);\n")
        f.write("\tmut_AAs = codonsToAminoAcids(mut_codons);\n")
        
        f.write("\t// get the mismatch AAs in the two sequences if there's any\n")
        f.write("\tif (any(mut_AAs!=original_AAs)) {\n")
        f.write("\t\t// if the site is conserved\n")
        f.write("\t\tif (match(mut.position, ConservSites)!=-1) {\n")
        f.write("\t\t\tmut.setSelectionCoeff(ConservCoeff);\n")
        f.write("\t\t}\n")
        
        f.write("\t\t//set selection coefficient based on buffering or not\n")
        f.write("\t\telse {\n")
        f.write("\t\t\tmut.setSelectionCoeff(0.0);\n")
        f.write("\t\t}\n")
        f.write("\t}\n")
        f.write("\treturn T;\n")
        f.write("}\n\n\n")
        
        
        f.write("// define a mutation callback during Wolbachia lost phase to control mutation generation:\n")
        f.write("// constrain the selection to similar amino acids closer to the ORIGINAL optimum\n")
        ## if wolbachia-infection phase comes first, then this section is defined to be the s2 section
        if (wolbachia_first):
            f.write("s2 %d:%d mutation(m1) {\n" % (neutral_burnIn_len+phase_len+1, neutral_burnIn_len+2*phase_len))
        ## if wolbachia-lost phase comes first, then this section is defined to be the s1 section
        else:
            f.write("s1 %d:%d mutation(m1) {\n" % (neutral_burnIn_len+1, neutral_burnIn_len+phase_len))
        f.write("\tsim.setValue(\"TotalMuts\", sim.getValue(\"TotalMuts\")+1);\n")
        f.write("\tsim_mut_aa_pos = integerDiv(sim.mutations.position, 3);\n")
        f.write("\tsim_mut_nuc_pos = c(sim_mut_aa_pos*3, sim_mut_aa_pos*3+1, sim_mut_aa_pos*3+2);\n")
        f.write("\t// if the mutation occurs at the codon already bearing a mutation and the site isn't conserved/neutral, block it\n")
        f.write("\tif (match(mut.position, sim_mut_nuc_pos)!=-1 & match(mut.position, ConservSites)==-1 & match(mut.position, NeutralSites)==-1){\n")
        f.write("\t\tsim.setValue(\"RejMuts\", sim.getValue(\"RejMuts\")+1);\n")
        f.write("\t\treturn F;\n")
        f.write("\t}\n")
        
        f.write("\t// if the mutation occurs at the neutral codon, return it with s = 0\n")
        f.write("\tif (match(mut.position, NeutralSites)!=-1){\n")
        f.write("\t\tmut.setSelectionCoeff(0.0);\n")
        f.write("\t\treturn T;\n")
        f.write("\t}\n")
        
        f.write("\t// get the original parental genome before mutation\n")
        f.write("\toriginal_genome = genome.nucleotides();\n")
        f.write("\t// get the new genome with mutation\n")
        f.write("\tmut_genome = substr(original_genome, 0, mut.position-1)+mut.nucleotide+substr(original_genome, mut.position+1);\n")
        f.write("\t// get the mutated and original codons\n")
        f.write("\toriginal_codons = nucleotidesToCodons(original_genome);\n")
        f.write("\tmut_codons = nucleotidesToCodons(mut_genome);\n")
        f.write("\toriginal_AAs = codonsToAminoAcids(original_codons);\n")
        f.write("\tmut_AAs = codonsToAminoAcids(mut_codons);\n")
        
        f.write("\t// get the mismatch AAs in the two sequences if there's any\n")
        f.write("\tif (any(mut_AAs!=original_AAs)) {\n")
        f.write("\t\t// if the site is conserved\n")
        f.write("\t\tif (match(mut.position, ConservSites)!=-1) {\n")
        f.write("\t\t\tmut.setSelectionCoeff(ConservCoeff);\n")
        f.write("\t\t}\n")
        
        f.write("\t\telse {\n")
        f.write("\t\t\tmismatch_idx = (mut_codons!=original_codons);\n")
        f.write("\t\t\toriginal_codon = original_codons[mismatch_idx];\n")
        f.write("\t\t\tmut_codon = mut_codons[mismatch_idx];\n")
        f.write("\t\t\tanc_codon = nucleotidesToCodons(AncSeq)[mismatch_idx];\n")
        f.write("\t\t\t// get the Miyata distance between the anc aa and the original aa and mut aa respectively\n")
        f.write("\t\t\tdiff_idx_mut = mut_codon*64 + anc_codon;\n")
        f.write("\t\t\tms_mut = M[diff_idx_mut];\n")
        f.write("\t\t\tdiff_idx_original = original_codon*64 + anc_codon;\n")
        f.write("\t\t\tms_original = M[diff_idx_original];\n")
        f.write("\t\t\t// if more diversifying compared with the original: deleterious\n")
        f.write("\t\t\tif (ms_mut - ms_original > 0)\n")
        f.write("\t\t\t\tmut.setSelectionCoeff(DelSelCoeff);\n")
        f.write("\t\t\t// if the same: neutral\n")
        f.write("\t\t\telse if (ms_mut - ms_original == 0)\n")
        f.write("\t\t\t\tmut.setSelectionCoeff(0.0);\n")
        f.write("\t\t\t// if more similar compared with the original (ms_mut-ms_original<0): beneficial\n")
        f.write("\t\t\telse\n")
        f.write("\t\t\t\tmut.setSelectionCoeff(BenSelCoeff);\n")
        f.write("\t\t}\n")
        f.write("\t}\n")
        f.write("\treturn T;\n")
        f.write("}\n\n\n")
        
        
        f.write("// create a population of 1000 dmel individuals\n")
        f.write("1 {\n")
        f.write("\tdefineConstant(\"simID\", getSeed());\n")
        f.write("\t// save this run's identifier, used to save and restore\n")
        f.write("\tc = sim.chromosome;\n")
        f.write("\tcatn(\"Ancestral: \" + c.ancestralNucleotides());\n")
        f.write("\tcatn('Ancestral Amino Acids: ' + codonsToAminoAcids(nucleotidesToCodons(c.ancestralNucleotides())));\n")
        f.write("\tdefineConstant(\"AncSeq\", sim.chromosome.ancestralNucleotides(format='char'));\n")
        f.write("\tsim.addSubpop('p1', 1000);\n\n")
        f.write("\t// Keep track the mutation generation process\n")
        f.write("\tsim.setValue(\"TotalMuts\", 0);\n")
        f.write("\tsim.setValue(\"RejMuts\", 0);\n")
        if (wolbachia_first):
            f.write("\t// register the wolbachia infection phases at the 1st and 3rd cycles\n")
        else:
            f.write("\t// register the wolbachia lost phases at the 1st and 3rd cycles\n")
        f.write("\tsim.rescheduleScriptBlock(s1, generations=c({}));\n".format(first_phases))
        if (wolbachia_first):
            f.write("\t// register the wolbachia lost phases at the 2nd and 4th cycles\n")
        else:
            f.write("\t// register the wolbachia infection phases at the 2nd and 4th cycles\n")
        f.write("\tsim.rescheduleScriptBlock(s2, generations=c({}));\n".format(second_phases))
        if (wolbachia_first):
            f.write("\t// register the entering wolbachia lost phases at the 2nd and 4th cycles\n")
        else:
            f.write("\t// register the entering wolbachia infection phases at the 2nd and 4th cycles\n")
        f.write("\tsim.rescheduleScriptBlock(s3, generations=seq({}, {}, {}));\n".format(neutral_burnIn_len+phase_len,wolbachia_non_burnIn_len+neutral_burnIn_len,phase_len*2))
        if (wolbachia_first):
            f.write("\t// register the entering wolbachia infection phases at the 1st and 3rd cycles\n")
        else:
            f.write("\t// register the entering wolbachia lost phases at the 1st and 3rd cycles\n")
        f.write("\tsim.rescheduleScriptBlock(s4, generations=seq({}, {}, {}));\n".format(neutral_burnIn_len+phase_len*2,wolbachia_non_burnIn_len+neutral_burnIn_len,phase_len*2))
        f.write("}\n\n\n")
        

        f.write("// Store the new reference sequence after neutral burn-in:\n")
        f.write("%d late() {\n" % neutral_burnIn_len)
        f.write("\tsim.chromosome.setAncestralNucleotides(AncSeq);\n")
        f.write("\tc = sim.chromosome;\n")
        f.write("\toriginal_genome = c.ancestralNucleotides();\n")
        f.write("\t// reset the new mutation's selection coefficient based on the Miyata distance it leads to the mutated amino acid\n")
        f.write("\tfor (mut in sim.mutations) {\n")
        f.write("\t\t// get the new genome with mutation\n")
        f.write("\t\tmut_genome = substr(original_genome, 0, mut.position-1)+mut.nucleotide+substr(original_genome, mut.position+1);\n")
        f.write("\t\t// get the mutated and original codons\n")
        f.write("\t\toriginal_codons = nucleotidesToCodons(original_genome);\n")
        f.write("\t\tmut_codons = nucleotidesToCodons(mut_genome);\n")
        f.write("\t\toriginal_AAs = codonsToAminoAcids(original_codons);\n")
        f.write("\t\tmut_AAs = codonsToAminoAcids(mut_codons);\n")
        
        f.write("\t\t// get the mismatch codon in the two sequences if there's any\n")        
        f.write("\t\t// Only if the mutation is nonsynonymous\n")
        f.write("\t\tif (any(mut_AAs!=original_AAs) & match(mut.position, ConservSites)==-1 & match(mut.position, NeutralSites)==-1) {\n")
        f.write("\t\t\tmismatch_idx = (mut_codons!=original_codons);\n")
        f.write("\t\t\tmut_codon = mut_codons[mismatch_idx];\n")
        f.write("\t\t\toriginal_codon = original_codons[mismatch_idx];\n")
        f.write("\t\t\t// get the Miyata distance between the translated aa's\n")
        f.write("\t\t\tdiff_index = mut_codon*64 + original_codon;\n")
        f.write("\t\t\tms_mut = M[diff_index];\n")
        if (wolbachia_first):
            f.write("\t\t\tmut.setSelectionCoeff(0.0);\n")
        else:
            f.write("\t\t\t// if more diversifying compared with the original: deleterious\n")
            f.write("\t\t\tif (ms_mut > 0)\n")
            f.write("\t\t\t\tmut.setSelectionCoeff(DelSelCoeff);\n")
            f.write("\t\t\t// if same: neutral\n")
            f.write("\t\t\telse\n")
            f.write("\t\t\t\tmut.setSelectionCoeff(0.0);\n")
        f.write("\t\t}\n")
        f.write("\t}\n")
        f.write("}\n\n\n")
        
        
        f.write("// entering the wolbachia lost phase: more neutral mutations lost 'protection' and become deleterious, more similar mutations to REF genome become beneficial\n")
        if (wolbachia_first):
            f.write("s3 %d late() {\n" % (neutral_burnIn_len+phase_len))
            f.write("\tcatn();\n")
            f.write("\tcatn(sim.generation+': Wolbachia Cycle Ends at '+(sim.generation-{})/{}+' AA Sequence: '+codonsToAminoAcids(nucleotidesToCodons(sim.chromosome.ancestralNucleotides())));\n".format(neutral_burnIn_len,phase_len))
        else:
            f.write("s4 %d late() {\n" % (neutral_burnIn_len+2*phase_len))
            f.write("\tcatn();\n")
            f.write("\tcatn(sim.generation+': No Wolbachia Cycle Ends at '+(sim.generation-{})/{}+' AA Sequence: '+codonsToAminoAcids(nucleotidesToCodons(sim.chromosome.ancestralNucleotides())));\n".format(neutral_burnIn_len,2*phase_len))
        
        f.write("\t// get the new reference genome after Wolbachia infection\n")
        f.write("\tref_genome = sim.chromosome.ancestralNucleotides();\n")
        f.write("\tref_codons = nucleotidesToCodons(ref_genome);\n")
        f.write("\tfor (mut in sim.mutations) {\n")
        f.write("\t\t// get the new genome with mutation\n")
        f.write("\t\tmut_genome = substr(ref_genome, 0, mut.position-1)+mut.nucleotide+substr(ref_genome, mut.position+1);\n")
        f.write("\t\t// get the mutated and original codons\n")
        f.write("\t\tmut_codons = nucleotidesToCodons(mut_genome);\n")
        f.write("\t\t// get the mismatch codon in the two sequences\n")
        f.write("\t\tif (any(mut_codons!=ref_codons) & match(mut.position, ConservSites)==-1 & match(mut.position, NeutralSites)==-1) {\n")
        f.write("\t\t\tmismatch_idx = (mut_codons!=ref_codons);\n")
        f.write("\t\t\tmut_codon = mut_codons[mismatch_idx];\n")
        f.write("\t\t\tref_codon = ref_codons[mismatch_idx];\n")
        f.write("\t\t\tanc_codon = nucleotidesToCodons(AncSeq)[mismatch_idx];\n")
        f.write("\t\t\t// get the Miyata distance between the translated aa's\n")
        f.write("\t\t\tdiff_index_mut = mut_codon*64 + anc_codon;\n")
        f.write("\t\t\tms_mut = M[diff_index_mut];\n")
        f.write("\t\t\tdiff_idx_ref = ref_codon*64 + anc_codon;\n")
        f.write("\t\t\tms_ref = M[diff_idx_ref];\n")
        f.write("\t\t\t// if more diversifying compared with the original: deleterious\n")
        f.write("\t\t\tif (ms_mut - ms_ref > 0)\n")
        f.write("\t\t\t\tmut.setSelectionCoeff(DelSelCoeff);\n")
        f.write("\t\t\t// if almost similar: neutral\n")
        f.write("\t\t\telse if (ms_mut - ms_ref == 0)\n")
        f.write("\t\t\t\tmut.setSelectionCoeff(0.0);\n")
        f.write("\t\t\t// if more similar compared with the original (ms_mut-ms_original<0): beneficial\n")
        f.write("\t\t\telse\n")
        f.write("\t\t\t\tmut.setSelectionCoeff(BenSelCoeff);\n")
        f.write("\t\t}\n")
        f.write("\t}\n")
        f.write("}\n\n\n")
        

        f.write("// enter the wolbachia infection phase: some previously deleterious mutations become neutral and similar mutations all become neutral\n")
        if (wolbachia_first):
            f.write("s4 %d late() {\n" % (neutral_burnIn_len+2*phase_len))
            f.write("\tcatn();\n")
            f.write("\tcatn(sim.generation+': No Wolbachia Cycle Ends at '+(sim.generation-{})/{}+' AA Sequence: '+codonsToAminoAcids(nucleotidesToCodons(sim.chromosome.ancestralNucleotides())));\n".format(neutral_burnIn_len,2*phase_len))
        else:
            f.write("s3 %d late() {\n" % (neutral_burnIn_len+phase_len))
            f.write("\tcatn();\n")
            f.write("\tcatn(sim.generation+': Wolbachia Cycle Ends at '+(sim.generation-{})/{}+' AA Sequence: '+codonsToAminoAcids(nucleotidesToCodons(sim.chromosome.ancestralNucleotides())));\n".format(neutral_burnIn_len,phase_len))
        
        f.write("\t// All mutations become neutral if not conserved\n")
        f.write("\tfor (mut in sim.mutations) {\n")
        f.write("\t\tif (match(mut.position, ConservSites)==-1 & match(mut.position, NeutralSites)==-1) {\n")
        f.write("\t\t\tmut.setSelectionCoeff(0.0);\n")
        f.write("\t\t}\n")
        f.write("\t}\n")
        f.write("}\n\n\n")
        
        
        f.write("// output the new ancestral sequence (with new substitutions) each 1k generations\n")
        f.write("// output the VCF file for 100 random samples each 1k generations\n")
        f.write("late() {\n")
        f.write("\tif ((sim.generation > 20000) & (sim.generation % 100 == 0) & ((sim.generation % 1000 != 0) & ((sim.generation-20000) % 3125 != 0))) {\n")
        f.write("\t\tcatn('Current mutations in generation '+sim.generation+':');\n")
        f.write("\t\tcat('{');\n")
        f.write("\t\tfor (mut in sim.mutations){\n")
        f.write("\t\t\tcat(mut.position+':'+mut.selectionCoeff+',');\n")
        f.write("\t\t}\n")
        f.write("\t\tcatn('}');\n")
        f.write("\t}\n\n")
        f.write("\tif ((sim.generation % 1000 == 0) | ((sim.generation-20000) % 3125 == 0)){\n")
        f.write("\t\tcatn();\n")
        f.write("\t\tcatn('Current mutations in generation '+sim.generation+':');\n")
        f.write("\t\tcat('{');\n")
        f.write("\t\tfor (mut in sim.mutations){\n")
        f.write("\t\t\tcat(mut.position+':'+mut.selectionCoeff+',');\n")
        f.write("\t\t}\n")
        f.write("\t\tcatn('}');\n")
        f.write("\t\t\n")
        f.write("\t\t// output to vcf\n")
        f.write("\t\tg = p1.sampleIndividuals({}).genomes;\n".format(pop_sample))
#         f.write("\t\tif ((sim.generation-20000) % 3125 == 0 & (sim.generation >= 20000)) {\n")
        f.write("\t\tif ((sim.generation-20000) % {} == 0 & (sim.generation >= 20000))\n".format(phase_len))
#         if (wolbachia_first):
        f.write("\t\t\tg.outputVCF('slim_bam_cds_buff_s{:.2e}_'+'Gen'+sim.generation+'.vcf', simplifyNucleotides=F);\n".format(sel_coeff))
#         else:
#             f.write("\t\t\tg.outputVCF('slim_bam_cds_buff_s'+BenSelCoeff+'_Gen'+sim.generation+'_wol-+_'+simID+'.vcf', simplifyNucleotides=F);\n")
        f.write("\t\tm = sortBy(unique(g.mutations), \"position\");\n")
        f.write("\t\tcatn('Generation '+sim.generation+\" segsites: \" + size(m));\n")
        f.write("\t\tpositions = m.position;\n")
        f.write("\t\tcatn('Generation '+sim.generation+\" positions: \" + paste(positions, \" \"));\n")
        f.write("\t\tnucleotides = m.nucleotide;\n")
        f.write("\t\tcatn('Generation '+sim.generation+\" nucleotides: \" + paste(nucleotides, \" \"));\n")
        f.write("\t\tcatn('Generation '+sim.generation+' Fixed mutations: ' + paste(sim.substitutions.nucleotide));\n")
        f.write("\t\tcatn('Their positions: ' + paste(sim.substitutions.position));\n")
        f.write("\t\t// print the new ancestral state with substitutions\n")
        f.write("\t\tc = sim.chromosome;\n")
        f.write("\t\tcatn('Generation '+sim.generation+' New Ancestral: ' + c.ancestralNucleotides());\n")
        f.write("\t\tcatn('Generation '+sim.generation+' AA Sequence: '+codonsToAminoAcids(nucleotidesToCodons(sim.chromosome.ancestralNucleotides())));\n")
        f.write("\t}\n")
        f.write("}\n\n\n")
        
        
        f.write("%d late() {\n" % (neutral_burnIn_len+wolbachia_non_burnIn_len))
        f.write("\t// print the fixed substitutions\n")
        f.write("\tcatn();\n")
        f.write("\tcatn('At the end of simulation in generation ' + sim.generation);\n")
        f.write("\tcatn('Fixed mutations: ' + paste(sim.substitutions.nucleotide));\n")
        f.write("\tcatn('Their positions: ' + paste(sim.substitutions.position));\n")
        f.write("\tcatn();\n\n")

        f.write("\t// print the new ancestral state with substitutions\n")
        f.write("\tc = sim.chromosome;\n")
        f.write("\tcatn('New Ancestral: ' + c.ancestralNucleotides());\n")
        f.write("\tcatn('New Ancestral Amino Acids: ' + codonsToAminoAcids(nucleotidesToCodons(c.ancestralNucleotides())));\n")
        f.write("\tcatn('Proportion of Rejected Mutations: ' + sim.getValue(\"RejMuts\")/sim.getValue(\"TotalMuts\"));\n")
        f.write("\tcatn();\n\n")
        f.write("\t// sim.outputFixedMutations();\n")
        f.write("}\n")
        

def get_opts_snpCalling(parser):
    ## input the ratio of neutral AAs
    parser.add_argument('--conservRatio', help='ratio of conserved AAs', required=True)
    ## input the ratio of selected AAs
    parser.add_argument('--selRatio', help='ratio of selected AAs', required=True)
    parser.add_argument('--sampleSize', help='sample size', required=True)
    ## input the selection coefficient
    parser.add_argument('--s', help='selection coefficient', required=True)
    parser.add_argument('--phaseLen', help='phase length', required=True)
    ## input the site type
    parser.add_argument('--run', help='sim run', required=True)
#     parser.add_argument('--sitesnum', help='total number of sites', required=True) ## deprecated
    return parser

        
def read_div_file(div_file):
    with open(div_file, 'r') as f:
        header = f.readline().split()
        subs_info = [int(x) for x in f.readline().split()]

    div_info = dict(zip(header, subs_info))
    return div_info


def read_mkt_2_pd(mkt_output_file):
    """
        mkt_output_file: the mkt output file from customized iMKT analysis
    """
    with open(mkt_output_file, 'r') as f:
        mkt_output = f.readlines()
    mkt_table = list(filter(lambda x: '$mktTable' in x, mkt_output))[0]
    mkt_idx = mkt_output.index(mkt_table)
    Ps, Ds = tuple(map(lambda x: int(x), mkt_output[mkt_idx:][2].strip().split()[2:4]))
    Pn, Dn = tuple(map(lambda x: int(x), mkt_output[mkt_idx:][3].strip().split()[2:4]))

    alpha_table = list(filter(lambda x: '$alphaCorrected' in x, mkt_output))[0]
    alpha_idx = mkt_output.index(alpha_table)
    cutoff0_alpha, cutoff0_pval = tuple(map(lambda x: float(x), mkt_output[alpha_idx:][2].strip().split()[2:]))
    cutoff05_alpha, cutoff05_pval = tuple(map(lambda x: float(x), mkt_output[alpha_idx:][3].strip().split()[2:]))
    cutoff10_alpha, cutoff10_pval = tuple(map(lambda x: float(x), mkt_output[alpha_idx:][4].strip().split()[2:]))

    return [Dn, Ds, Pn, Ps, cutoff0_alpha, cutoff0_pval, cutoff05_alpha, cutoff05_pval, cutoff10_alpha, cutoff10_pval]


def main():
    buff_prefix = "slim_bam_cds_buff_"
    start_time = time.time()
    
    conserv_coeff = -0.1
    neutral_burnIn_len = 20000
    wolbachia_non_burnIn_len = 25000
#     wolbachia_first = bool(int(sys.argv[3]))
    wolbachia_first = True
    print("Wolbachia first:", wolbachia_first)
    
    # Initiate the parser
    parser = argparse.ArgumentParser()
    # Add long and short argument
    parser = get_opts_snpCalling(parser)
    # Read arguments from the command line
    args = parser.parse_args()
    
    sel_coeff = float(args.s)
    phase_len = int(args.phaseLen)
    run_number = int(args.run)
    conserv_aa_ratio = float(args.conservRatio)
    sel_aa_ratio = float(args.selRatio)
    pop_sample = int(args.sampleSize)
    
    print('SLiM Buff Base Model with ConservRatio={} SelRatio={} S={} s={} phase={} Run{}'.format(conserv_aa_ratio, sel_aa_ratio, pop_sample, sel_coeff, phase_len, run_number))
    
    wolbachia_seq = '+-' if wolbachia_first else '-+'
    buff_file = buff_prefix + 's{:.0e}_phase{}k_wol{}_run{}'.format(sel_coeff, str(phase_len/1000).replace('.','pt'), wolbachia_seq, run_number)
    
    wrt2buff(buff_file+'.slim', conserv_aa_ratio, sel_aa_ratio, conserv_coeff, pop_sample, sel_coeff, phase_len, neutral_burnIn_len, wolbachia_non_burnIn_len, wolbachia_first)
    
    print("Successfully output files at {}".format(buff_file))
    print("Run", buff_file)
    
    with open(buff_file+'.out', "w+") as output:
        subprocess.call(["slim", buff_file+'.slim'], stdout=output);
        
    print("SLiM Run for {}s".format(time.time()-start_time))
    
    vcf_file = 'slim_bam_cds_buff_s{:.2e}_Gen{}.vcf'.format(sel_coeff, neutral_burnIn_len+wolbachia_non_burnIn_len)
    fasta_file = '.'.join(vcf_file.split('.')[:-1])+'.fa'
    daf_file, div_file = '.'.join(vcf_file.split('.')[:-1])+'_2fold.daf', '.'.join(vcf_file.split('.')[:-1])+'_2fold.div'
    mkt_output_file = buff_file + '_imkt.out'
    
    with open(buff_file+'.out', 'r') as f:
        f.readline()
        slim_seed = f.readline().strip()
    
    with open('vcf2fasta.out', 'w') as output:
        subprocess.run(["python", "vcf2fasta.py", vcf_file, buff_file+'.out'], stdout=output)
        
    with open('sfsFromFasta_2fold.out', 'w') as output:
        subprocess.run(["python", "sfsFromFasta_2fold.py", 
                         fasta_file, daf_file, div_file], stdout=output)
        
    div_info = read_div_file(div_file)
    print(div_info)
    
    with open(mkt_output_file, 'w') as output:
        subprocess.run(["/programs/R-4.0.0/bin/Rscript", "--vanilla",
                         "run_iMKT.R", daf_file, div_file], stdout=output)
        
    imkt_results = read_mkt_2_pd(mkt_output_file)
    
    simulation_runs = '{}_params.txt'.format(buff_file)
    with open(simulation_runs, 'a') as f:
        f.write('Run{}:\nConservRatio={}\tSelRatio={}\tS={}\ts={}\tseed={}\tDi={}\tD0={}\tmi={}\tm0={}\n'.format(run_number, conserv_aa_ratio, sel_aa_ratio, pop_sample, sel_coeff, slim_seed, div_info['Di'], div_info['D0'], div_info['mi'], div_info['m0']))
        f.write('Dn={}\tDs={}\tPn={}\tPs={}\talpha0={}\tpval0={}\talpha05={}\tpval05={}\talpha10={}\tpval10={}\n\n'.format(*imkt_results))
    
    
if __name__ == '__main__':
    main()