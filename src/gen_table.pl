#!/usr/bin/perl -w

@genome_list=("A1", "A2", "A3");
@error_list=("0p0", "0p2", "0p5", "2p3");
@length_list=("075","100","150","200","300","400","600","1000");

for $genome (@genome_list) {
    for $error (@error_list) {
        for $length (@length_list) {
            $stem = $genome.'-'.$error.'-'.$length;
            $N = `confusion.pl $stem`;
            print $N;
        }
    }
}