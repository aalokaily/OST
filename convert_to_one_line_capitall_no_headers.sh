
for file in $(find /gcrf/KHCC/CTAG/NGS/Compression/paper_test_dataset/fasta_files -maxdepth 2 -name "*.fasta" ); do
file_name=$(basename $file)
file_name=${file_name:0:-6}

echo $file_name

python /gcrf/KHCC/CTAG/NGS/scripts/filter_DNA_file_to_4_bases_and_N.py $file
#grep -v '^>' $file | tr -d '"\r\nN' | tr '[a-z]' '[A-Z]' > $file_name.one_line
#grep -v '^>' $file | tr -d '"\r\n' | tr '[a-z]' '[A-Z]' > $file_name.one_line

done

