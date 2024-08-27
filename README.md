# OST 
The scripts are implementation of OST algorithm which is a novel encoding algorithm for data compression. The full description is available at https://www.biorxiv.org/content/10.1101/2020.08.24.264366v1.full. The name of encoding algorithm is after Anas Al-okaily, Pramod Srivastava, and Abdelghani Tbakhi (Okaily-Srivastava-Tbakhi).

Briefly and generally speaking, the algorithm starts by scanning the DNA data with non-overlapping windows, label the sequence within each window, concatenate it with the sequences in a bin correspondent to that label, and then output the label into a stream. Now, encode the labels of the bins based on, for instance, the number of the sequences in the bins; compress the stream of labels using the label codes. Then, compress the sequences in each bin using suitable compression algorithm depending on the content of the sequences in that bin. Finally, compress all the compression results (bins and stream of labels) together. 

The decompression process will be by firstly decompressing each bin and the stream of labels. Then, read the labels sequentially from the stream of labels, at each reading and using a counter for each bin, get the next sequence (of length same as the length of the non-overlapping window used during the compression process) from the bin of that label then increment the counter at that bin.

The seven scripts OST-DNA-bcm, OST-DNA-brotli, OST-DNA-bsc, OST-DNA-bzip2, OST-DNA-lrzip, OST-DNA-lzip, and OST-DNA-xz are python scripts for applying OST algorithm on DNA data and where the correspondent tools (bcm, brotli, bsc, bzip2, lrzip, lzip, and xz) are used to compress the bins. For each tool a window length and label length are parameters. The description of these parameters are available at https://www.biorxiv.org/content/10.1101/2020.08.24.264366v1.full .

For each script and in order to validate the command used for compression the bins (line number 149 in each script) and decompression the bins (line number 233 in each script) or change them according to your need, please go inside the script to do so.

----------------------------------------------- Preparation --------------------------------------------------------------------------------------------

Firstly, you may convert the genome in fasta format to a one-line genome which remove any non A, C, G, T, and N (case is sensitive) and headers. This can be done using the script filter_DNA_file_to_4_bases_and_N.py by running the command:

python filter_DNA_file_to_4_bases_and_N.py $file.fasta 

----------------------------------------------- Running the scripts ------------------------------------------------------------------------------------
Compression

python OST-DNA-xxx.py $genome $label_length $window_length 

Note: huffman package must be installed, as the scripts import this package. 

Decompression 

python OST-DNA-xxx.py -d $genome.ost.7z


For contact, please email AA.12682@KHCC.JO (the email of the first author Anas Al-okaily).
