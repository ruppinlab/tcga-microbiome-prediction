#!/usr/bin/env perl

use strict;
use warnings;
use autodie;
use FindBin;
use File::Basename qw(fileparse);
use File::Path qw(make_path remove_tree);
use List::Util qw(any);

chdir $FindBin::Bin;

my $data_dir = "$FindBin::Bin/data";

# survival
my @surv_files = glob("$data_dir/tcga_*_surv_*_eset.rds");
for my $file (@surv_files) {
    my $file_basename = fileparse($file, qr/\.[^.]*/);
    my @file_basename_parts = split('_', $file_basename);
    my $cancer = $file_basename_parts[1];
    my $cancer_target = join('_', @file_basename_parts[1, 3]);
    my $data_type = $file_basename_parts[4];

    my $lscratch = 400;
    my $mem_per_cpu = '3584m';
    my $num_jobs = $data_type eq 'kraken' ? 32 : 56;
    my $scv_splits = 4;
    my $scv_repeats = 5;
    my $test_splits = 100;
    my $test_size = 0.25;

    if (any { $cancer_target eq $_ } qw(
        chol_os  chol_pfi
                 dlbc_pfi
        kich_os
        thym_os
    )) {
        $scv_splits = 3;
        $scv_repeats = 5;
    } elsif (any { $cancer_target eq $_ } qw(
        dlbc_os
        pcpg_os
        tgct_os
    )) {
        $scv_splits = 2;
        $scv_repeats = 5;
    }

    my $results_dir = 'results/' . join(
        '_', @file_basename_parts[0 .. $#file_basename_parts - 1], 'cnet'
    );
    remove_tree($results_dir, { safe => 1 }) if -d $results_dir;
    make_path($results_dir);

    my $cmd = <<"    CMD";
    ./submit_slurm_model.sh
    --dataset $file
    --scv-splits $scv_splits
    --scv-repeats $scv_repeats
    --test-splits $test_splits
    --test-size $test_size
    --n-jobs $num_jobs
    --verbose 1
    --scv-verbose 0
    --out-dir $results_dir
    --sbatch-opts "--gres=lscratch:$lscratch --mem-per-cpu=$mem_per_cpu --time=4-00:00:00 --output=$results_dir/slurm-\%j.out"
    CMD
    $cmd =~ s/^\s+//;
    $cmd =~ s/\s+$//;
    $cmd =~ s/\s+/ /g;
    system($cmd);
}

# response
my @kraken_resp_files = glob("$data_dir/tcga_*_resp_*_eset.rds");
for my $file (@kraken_resp_files) {
    my $file_basename = fileparse($file, qr/\.[^.]*/);
    my @file_basename_parts = split('_', $file_basename);
    my $cancer = $file_basename_parts[1];
    my $cancer_target = join('_', @file_basename_parts[1, 3]);
    my $data_type = $file_basename_parts[4];

    my $lscratch = 400;
    my $mem_per_cpu = '3584m';
    my $num_jobs = $data_type eq 'kraken' ? 32 : 56;
    my $scv_splits = 3;
    my $scv_repeats = 5;
    my $test_splits = 4;
    my $test_repeats = 25;

    if (any { $cancer_target eq $_ } qw(
        stad_oxaliplatin
    )) {
        $test_splits = 3;
        $test_repeats = 33;
    }

    my $results_dir = "results/$file_basename_parts[2]/" . join(
        '_', @file_basename_parts[0 .. $#file_basename_parts - 1], 'rfe'
    );
    remove_tree($results_dir, { safe => 1 }) if -d $results_dir;
    make_path($results_dir);

    my $cmd = <<"    CMD";
    ./submit_slurm_model.sh
    --dataset $file
    --scv-splits $scv_splits
    --scv-repeats $scv_repeats
    --test-splits $test_splits
    --test-repeats $test_repeats
    --n-jobs $num_jobs
    --verbose 1
    --scv-verbose 0
    --out-dir $results_dir
    --sbatch-opts "--gres=lscratch:$lscratch --mem-per-cpu=$mem_per_cpu --time=4-00:00:00 --output=$results_dir/slurm-\%j.out"
    CMD
    $cmd =~ s/^\s+//;
    $cmd =~ s/\s+$//;
    $cmd =~ s/\s+/ /g;
    system($cmd);
}
