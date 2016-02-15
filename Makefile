# ***********************************************************************
#
#  Makefile for SingleSplice
#  =====================================
#  Josh Welch
#  2/13/16
# ***********************************************************************
SAMTOOLS=./diffsplice/src/parse_bam/samtools-0.1.18/
BOOST=./boost_1_60_0/
.PHONY: clean

all: test_ratio_change samtools parse_bam parse_frag diff_gtree diff_asm_analysis diff_expr_analysis diffsplice_exec

test_ratio_change: test_ratio_change.o
	g++ test_ratio_change.o -o test_ratio_change
test_ratio_change.o: test_ratio_change.cc
	g++ -c -O3 test_ratio_change.cc -o test_ratio_change.o -I${BOOST}
samtools: 
	${MAKE} -C ${SAMTOOLS}
parse_bam: diffsplice/src/parse_bam/parse_bam.cpp
	g++ -I${SAMTOOLS} -g -O2 -Wall ./diffsplice/src/parse_bam/parse_bam.cpp -o ./diffsplice/bin/parse_bam -lz -L${SAMTOOLS} -lbam
parse_frag: diffsplice/src/parse_frag/parse_frag.cpp
	g++ -m64 -o ./diffsplice/bin/parse_frag ./diffsplice/src/parse_frag/parse_frag.cpp
diff_gtree: diffsplice/src/diff_gtree_ref/common_function.cpp diffsplice/src/diff_gtree_ref/decomposition.cpp diffsplice/src/diff_gtree_ref/estimation.cpp diffsplice/src/diff_gtree_ref/input_data.cpp diffsplice/src/diff_gtree_ref/output.cpp diffsplice/src/diff_gtree_ref/splice_graph.cpp diffsplice/src/diff_gtree_ref/cut_exon_bound.cpp diffsplice/src/diff_gtree_ref/wavelet2s.cpp diffsplice/src/diff_gtree_ref/main.cpp
	g++ -m64 -I ./diffsplice/src/diff_gtree_ref/wavelet/include/ -o ./diffsplice/bin/diff_gtree diffsplice/src/diff_gtree_ref/common_function.cpp diffsplice/src/diff_gtree_ref/decomposition.cpp diffsplice/src/diff_gtree_ref/estimation.cpp diffsplice/src/diff_gtree_ref/input_data.cpp diffsplice/src/diff_gtree_ref/output.cpp diffsplice/src/diff_gtree_ref/splice_graph.cpp diffsplice/src/diff_gtree_ref/cut_exon_bound.cpp diffsplice/src/diff_gtree_ref/wavelet2s.cpp diffsplice/src/diff_gtree_ref/main.cpp -L diffsplice/src/diff_gtree_ref/wavelet/lib/ -lfftw3 -lm -std=c++0x 
diff_asm_analysis: diffsplice/src/diff_asm_analysis/asm_def.cpp diffsplice/src/diff_asm_analysis/calc_function.cpp diffsplice/src/diff_asm_analysis/difftrans_test.cpp diffsplice/src/diff_asm_analysis/input_data.cpp diffsplice/src/diff_asm_analysis/permutation_test.cpp diffsplice/src/diff_asm_analysis/pairwise_test.cpp diffsplice/src/diff_asm_analysis/write_output.cpp diffsplice/src/diff_asm_analysis/main.cpp
	g++ -m64 -o ./diffsplice/bin/diff_asm_analysis ./diffsplice/src/diff_asm_analysis/asm_def.cpp ./diffsplice/src/diff_asm_analysis/calc_function.cpp ./diffsplice/src/diff_asm_analysis/difftrans_test.cpp ./diffsplice/src/diff_asm_analysis/input_data.cpp ./diffsplice/src/diff_asm_analysis/permutation_test.cpp ./diffsplice/src/diff_asm_analysis/pairwise_test.cpp ./diffsplice/src/diff_asm_analysis/write_output.cpp ./diffsplice/src/diff_asm_analysis/main.cpp
diff_expr_analysis: ./diffsplice/src/diff_expr_analysis/gene_def.cpp diffsplice/src/diff_expr_analysis/calc_function.cpp diffsplice/src/diff_expr_analysis/difftrans_test.cpp diffsplice/src/diff_expr_analysis/input_data.cpp diffsplice/src/diff_expr_analysis/permutation_test.cpp diffsplice/src/diff_expr_analysis/pairwise_test.cpp diffsplice/src/diff_expr_analysis/write_output.cpp diffsplice/src/diff_expr_analysis/main.cpp
	g++ -m64 -o ./diffsplice/bin/diff_expr_analysis ./diffsplice/src/diff_expr_analysis/gene_def.cpp ./diffsplice/src/diff_expr_analysis/calc_function.cpp ./diffsplice/src/diff_expr_analysis/difftrans_test.cpp ./diffsplice/src/diff_expr_analysis/input_data.cpp ./diffsplice/src/diff_expr_analysis/permutation_test.cpp ./diffsplice/src/diff_expr_analysis/pairwise_test.cpp ./diffsplice/src/diff_expr_analysis/write_output.cpp ./diffsplice/src/diff_expr_analysis/main.cpp
diffsplice_exec: ./diffsplice/src/diffsplice/diffsplice.cpp
	g++ -m64 -o ./diffsplice/bin/diffsplice ./diffsplice/src/diffsplice/diffsplice.cpp

clean:
	${MAKE} -C ${SAMTOOLS} clean
	rm test_ratio_change test_ratio_change.o ./diffsplice/bin/*
