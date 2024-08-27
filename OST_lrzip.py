import sys 
import os
import math
import random
import fileinput
import huffman
import subprocess
from datetime import datetime
from bitstring import BitArray
from collections import defaultdict
from collections import OrderedDict 
import operator
import ast
import time



def compute_windows_kmers(original_file, window_length, k):
	temp_dict = defaultdict(defaultdict)    	
	temp_dict.clear()
    	files_dict = defaultdict()
    	tree_keys_file = open(os.path.basename(original_file).rsplit(".",1)[0] + "." + os.path.basename(original_file).rsplit(".",1)[1] + ".ost.txt", 'w')
	with open(original_file) as f:
            while True:
                window = f.read(window_length)
                if not window:
                    break
                else:
                    current_length = len(window)  # in case l is not multiple of current_length
                    if current_length > 0:
                        #collect freq of each char to classify the window and name the tree
                        chars_freq_list = []
                        chars_freq_dict = defaultdict(int)
                        chars_freq_dict.clear()
                        for char in window:
                            if char in chars_freq_dict:
                                chars_freq_dict[char] += 1
                            else:
                                chars_freq_dict[char] = 1

                    # Add the char and their freqencies to a list format so that it can be inputted to Huffman encoding function 
                    for char in chars_freq_dict:
                        chars_freq_list.append((char, chars_freq_dict[char]))
                    
                    # build the key from the huffman tree 
                    tree_key = ""
                    chars_freq_list = sorted(chars_freq_list, key = lambda x: x[1], reverse=True)
                    t = huffman.codebook(chars_freq_list)
                    for key in sorted(t.items(), key=lambda x: len(x[1]))[:label_length]:
                        tree_key += key[0]
                    tree_key += "_"
                    for code in sorted(t.items(), key=lambda x: len(x[1]))[:label_length]:
                        if len(code[1]) == 0:  # this due to the bug in huffman package that if the input is one character then this package assign no code for this only character
                            tree_key += str(1)
                        else:
                            tree_key += str(len(code[1]))

                    # write seqeunces of same tree key to a seperate file with same name of tree key
                    if tree_key in files_dict:
                        files_dict[tree_key].write(window)
                    else:
                        files_dict[tree_key] = open(tree_key+".txt", "w")
                        files_dict[tree_key].write(window)                    

                    
                    # compute the freqeuncies of the tree keys
                    if tree_key in temp_dict["tree_key_freq"]:
                        temp_dict["tree_key_freq"][tree_key] += 1
                    else:
                        temp_dict["tree_key_freq"][tree_key] = 1
                    
                    # write the tree_key to a file stream
                    tree_keys_file.write(tree_key + "\n")
			
	return temp_dict

def generate_huffman_trees_for_windows_kmers(all_huffmans_trees):
	# compute the saving that is made by each tree in order to use this saving in
	# building the huffman tree for tree keys so that the tree key with larger saving
	# get smaller codeword
	temp = [(u, v) for u, v in all_huffmans_trees["tree_key_freq"].items()]
	huffman_for_trees = huffman.codebook(temp)

    	if len(temp) == 1:
        	all_huffmans_trees["tree_code"][temp[0][0]] = "0"
    	else:
        	for key in sorted(all_huffmans_trees["tree_key_freq"].items(), key=lambda x: x[1], reverse=True):
            		all_huffmans_trees["tree_code"][key[0]] = huffman_for_trees[key[0]]

	return all_huffmans_trees

def compress_DNA(original_file, final_huffmans_trees, window_length, start_time):
    # convert and store the binary string to binary original_file
    # add needed data at the begining of binary original_file
    binary_file = open(os.path.basename(original_file).rsplit(".",1)[0] + "." + os.path.basename(original_file).rsplit(".",1)[1] + ".ost", 'w')

    out_string = ""
    binary_file.write("Left_over_bits::########" + "\n")
    binary_file.write("Window_length::" + str(window_length))
    binary_file.write("\nLabel_length::" + str(label_length))
    binary_file.write("\ntree_code_dictionary::" + str(dict(final_huffmans_trees["tree_code"])))
    #binary_file.write("\nkmers_codewords_dictionary::" + str(dict(final_huffmans_trees["kmers_codewords"])))
    binary_file.write("\n#####################") 
    
    ll = 0
    tree_keys_file = open(os.path.basename(original_file).rsplit(".",1)[0] + "." + os.path.basename(original_file).rsplit(".",1)[1] + ".ost.txt", "rw+")
    while True:
        tree_key = tree_keys_file.readline().strip()
        #print tree_key
        if not tree_key:
            break
        else:  
            tree_code = final_huffmans_trees["tree_code"][tree_key]
            out_string += tree_code
            
            ii = 0
            ll += len(tree_code)
            if ll > 8000:
                while ll >= 8:
                    binary_file.write(chr(int(out_string[ii:ii+8], 2)))
                    #print out_string[ii:ii+8], chr(int(out_string[ii:ii+8], 2))
                    ll = ll - 8
                    ii += 8
                out_string = out_string[ii:]

    #encode the rest of the data
    ii = 0
    ll = len(out_string)
    while ll >= 8:
        binary_file.write(chr(int(out_string[ii:ii+8], 2)))
        ll = ll - 8
        ii += 8
    out_string = out_string[ii:]
    
    binary_file.seek(0,0)
    binary_file.write("Left_over_bits::" + out_string)
    binary_file.close() 
    os.system("rm " + os.path.basename(original_file).rsplit(".",1)[0] + "." + os.path.basename(original_file).rsplit(".",1)[1] + ".ost.txt")
    
    # compress the bins
    original_file_size = int(os.path.getsize(original_file))
    original_file_name = os.path.basename(original_file).rsplit(".",1)[0] + "." + os.path.basename(original_file).rsplit(".",1)[1]
    s = 0
    tree_keys_by_size = sorted(final_huffmans_trees["tree_key_freq"].items(), key=lambda x: x[1], reverse=True)
    temp = [i[0] for i in tree_keys_by_size]  # to compress larger file first
    
    os.system("echo start logging > " + original_file_name + ".log")
    for tree_key in temp:
       	os.system("lrzip -q -N 1 -f " + tree_key + ".txt -o " + tree_key + ".txt.lrz &>> " + original_file_name + ".log")
        print "window_length", window_length, "label_length", label_length, tree_key, os.path.getsize(tree_key + ".txt"), os.path.getsize(tree_key + ".txt.lrz"), float(os.path.getsize(tree_key + ".txt.lrz"))/os.path.getsize(tree_key + ".txt")
       	os.system("rm " + tree_key + ".txt")
        os.system("7za a " + original_file_name + ".ost.7z " + tree_key + ".txt.lrz &>> " + original_file_name + ".log")
      	s += os.path.getsize(tree_key + ".txt.lrz")
        os.system("rm " + tree_key + ".txt.lrz")

    
    s += os.path.getsize(original_file_name + ".ost")
    print "window_length", window_length, "label_length", label_length, "without 7z: CR ", round(float(s)/original_file_size, 4), "CT", int(time.time() - start_time), " seconds"
    
    os.system("7za a " + original_file_name + ".ost.7z " + original_file_name + ".ost &>> " + original_file_name + ".log")
    d = float(os.path.getsize(original_file_name + ".ost.7z"))/original_file_size
    print "window_length", window_length, "label_length", label_length, "with    7z: CR ", round(d, 4), "CT", int(time.time() - start_time), " seconds"
    
    os.system("rm " + original_file_name + ".ost")
    #os.system("rm " + original_file_name + ".log")
    return True

def compress(file, window_length, label_length):    
    start_time = time.time() 
    all_huffmans_trees = compute_windows_kmers(file, window_length, label_length)
    final_huffmans_trees = generate_huffman_trees_for_windows_kmers(all_huffmans_trees)
    compress_DNA(file, final_huffmans_trees, window_length, start_time)


    
    
    
    
###########################   Main   ####################################################################

if sys.argv[1] != "-d":

    file = sys.argv[1]    
    label_length = int(sys.argv[2])    
    window_length = int(sys.argv[3])

    compress(file, window_length, label_length)
    

############ Decompression ###########################      
else:
    if len(sys.argv) != 3:
        print "please provide OST file to uncompress.\n"
        
    OST_file = sys.argv[2]
    
    start_time = time.time() 
    os.system("7za e " + OST_file + " &> log.log")     
    labels_stream_bin_file = subprocess.check_output("ls *.ost", shell=True).strip()
    os.system("mv log.log " + labels_stream_bin_file + ".log")
    
    with open(labels_stream_bin_file) as f:
        header = ""
        while True:
            window = f.read(1)
            if not window:
                break
            else:
                header += window

            if header[-21:] == "#####################":
                break
   
    for line in header.split("\n"):
        if line.split("::")[0] == "Window_length":
            Window_length = int(line.split("::")[1])
            
        if line.split("::")[0] == "Label_length":
            label_length = int(line.split("::")[1])

        if line.split("::")[0] == "tree_code_dictionary":
            tree_code_dictionary = ast.literal_eval(line.split("::")[1])
            # swap key, value in dictionary
            tree_code_dictionary = {value:key for key, value in tree_code_dictionary.items()}

        if line.split("::")[0] == "Left_over_bits":
            Left_over_bits = line.split("::")[1].split("#")[0].strip() 
    
    # start the decompression process
    # decompress the bins
    for key in tree_code_dictionary:
        #print tree_code_dictionary[key]
        os.system("lrzip -d -q -N 1 -f " + tree_code_dictionary[key] + ".txt.lrz -o " + tree_code_dictionary[key] + ".txt &>> " + labels_stream_bin_file + ".log")      
        os.system("rm " + tree_code_dictionary[key] + ".txt.lrz" + " &>> " + labels_stream_bin_file + ".log")
        globals()[tree_code_dictionary[key]] = open(tree_code_dictionary[key] + ".txt", 'r')
    
    
    # decompress the labels_stream
    decompressed_file = open(labels_stream_bin_file[:-4], 'w')
    with open(labels_stream_bin_file) as f:
            #print header
            f.read(len(header))
            bits_string = ""
            last_bits = 0
            while True:
                chars = f.read(1000)
                if not chars:
                    bits_string += Left_over_bits
                    last_bits = 1                    
                else:
                    for c in chars:
                        bits_string += bin(ord(c))[2:].zfill(8)
                        #print bin(ord(c))[2:].zfill(8), c
                    if len(chars) < 1000:
                        bits_string += Left_over_bits
                        last_bits = 1 
                        
                l = len(bits_string)                
                i = 0                    
                while i < l:
                    if i > 8 * 900 and last_bits == 0:
                        bits_string = bits_string[i:]
                        break
                        
                    tree_code_decoded = True
                    j = 1
                    while tree_code_decoded and j <= l:
                        tmp = str(bits_string[i:i+j])
                        #print tmp, j, i, l
                        if tmp in tree_code_dictionary:
                            tree_code = tree_code_dictionary[tmp]
                            decompressed_seq = globals()[tree_code].read(Window_length)
                            #print tree_code, tmp, j, i, l, decompressed_seq
                            decompressed_file.write(decompressed_seq)
                            tree_code_decoded = False                            
                        else:
                            j += 1                    
                    i += j
                    
                if last_bits == 1:
                    break
    
    decompressed_file.close()
    
    for key in tree_code_dictionary:
        os.system("rm " + tree_code_dictionary[key] + ".txt")
    
    os.system("rm " + labels_stream_bin_file)
    print "window_length", Window_length, "label_length", label_length, "Decompression time is ", int(time.time() - start_time), " seconds"

        



