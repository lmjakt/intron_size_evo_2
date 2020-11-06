#!/bin/bash

find $HOME/genomes/introns/genome_sequences_98/ -name "*.nsq" | sed -E 's/\.?[0-9]*\.nsq$//' | sort -u
