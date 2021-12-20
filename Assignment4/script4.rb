require 'bio'

#This script is think only for search orthologues between two species
#The usaege is: ruby script4.rb file1 (Spombe.fa) file2 (Tair.fa)
#The script creates a report file in which the ID of the orthologues, the script runtime and the number of orthologues is provided

#Method for calculating the running time of the script 
#This method is implemented from https://stackoverflow.com/questions/11406410/measure-and-benchmark-time-for-ruby-methods
@starting = Process.clock_gettime(Process::CLOCK_MONOTONIC)

#Importing the proteome files of the species which are going to be compared to obtain the orthologues between them
#First, I recommend to introduce the S. Pombe file as its proteome has less sequences than Tair's, so the time of execution will be reduced
file1 = ARGV[0]
file2 = ARGV[1]


#This function creates two different lists, one for each proteome file of input. This list include the input file name (which es the same for the database created for the blast),
#the type of the sequences that appear in the input file; in this case, the type can be nucleotide ('nucl') or protein ('prot'), and
#the flat file created with the bioruby method [Bio::FlatFile.auto()] for each input file

def important_info (file1 = file1, file2 = file2)
	#Create the two arrays
	@file1_list = Array.new
	@file2_list = Array.new

	#Create a list ('files') of the two files of input (proteome files) to iterate over it and obtain the wanted outputs
	files = [file1, file2]
	#Counter to not overwrite the lists in the loop
	contador = 0
	#loop to iterrate over the list 'files'
	files.each do |file|
		#Obtain the flat file from the proteome files
		file_flat = Bio::FlatFile.auto(file)
		#Obtain the sequence types of the porteome files
		file_type = file_flat.next_entry.to_biosequence.guess
		if file_type == Bio::Sequence::NA
			f_type_ = 'nucl'
		elsif file_type == Bio::Sequence::AA
			f_type_ = 'prot' 
		end
		
		#Create the command to obtain the databases that are needed to do the reciprocal blasts, this command needs an input fasta file (proteomes files)
		#and the type of sequences in the file, it can be provided also a name for the output file (the database) but this is optative
		#if there is no name given to blast, the new file is titled as the input proteome file
		system("makeblastdb -in #{file} -dbtype '#{f_type_}'")
	
		contador +=1
		#Store the outputs in the lists for each input file
		if contador == 1
			@file1_list = [file,f_type_,file_flat]

		else
			@file2_list = [file,f_type_,file_flat]
	
		end
	end
	#The function returns both list which have the information that would be necessary in the next steps
	return @file1_list
	return @file2_list

end

#Calling the function
important_info(file1,file2)

#Function to search blast type
#This function needs the two arrays created by the first funtion as input
def search_blast_type(data1 = @file1_list, data2 = @file2_list)
	#Get the type of sequences from each file from the arrays, this feature is in the second [1] possition of the array
	type1 = data1[1]
	type2 = data2[1]

	#if the query sequence is a nucleotide sequence and the database against which the blast is performed is of aa sequences --> blastx
	if type1 == "nucl" && type2 == "prot"
		blastype = "blastx"
	#if the query sequence is a protein sequence and the database against which the blast is performed is of nucleotide sequences --> tblastn
	elsif type1 == "prot" && type2 == "nucl"
		blastype = "tblastn"
	#if the query sequence is a protein sequence and the database against which the blast is performed is of aa sequeneces  --> blastp
	elsif type1 == "prot" && type2 == "prot"
		blastype = "blastp"
	#if the query sequence is a nucleotide sequence and the database against which the blast is performed is of nucleotide sequences --> blastn
	elsif type1 == "nucl" && type2 == "nucl"
		blastype = "blastn"
	end
	#The function returns the type of blast
	return blastype
end

#Function to perform blast
#To run a blast you need a query sequence, a database and the type of blast
#the function returns the sequence ID of the best hit for each query sequence (that is porvided as input) against the database provided 
def blast(blast_type = nil, query = nil , db = nil)
	#This is the coverage filter to consider that the blast result is good enough to be taken into consideration
	coverage_limit = 50
	#This step is the blast, in here we need the blast type, the database and filtering (-F T -s F)  and e-value arguments (-e 1e-6)
	factory = Bio::Blast.local("#{blast_type}", "#{db}","-F T -s F  -e 1e-6") 
	#Obtaining the blast report to only select the first result (the best)
	report = factory.query(query)
	#if the report is empty
	if report.hits[0] == nil 
		#This is for me to know in which step of the analysis is the script, so it would appear in the screen if there is no resut in the blast. 
		#The type, the query secuence and the database that is used in the blast is printed in the screen.
		puts "Report is empty for the query secuence " + "#{query.entry_id}" + " when running a " + "#{blast_type}" + " against " + "#{db.split("/")[5]}" + " database."
	else
		#if there are results
		seq = report.hits[0]
		#COVERAGE ANALYSIS
		
		#coverage = [(end positions of the query sequence in the alignment) - (start positions of the query sequence in the alignment)] / lenth of the query sequence 
 		if blast_type == "blastx"
 			#if query sequence is a nucleotide sequence and the database against which the blast is performed is of aa sequences, that means we are performing a blastx,
			#the length of the query sequence has to be divided between 3 as 3 nucleotides is an amino acid)
 			#calculate the coverage of the hit sequence and and multiply it by 100 as the coverage limit is of 50 not 0.5
 			coverage_s = ((seq.query_end.to_f - seq.query_start.to_f)/(seq.query_len.to_f/3)) * 100

 		else 
 		#if we are not in blastx case the calculation of coverage is simple:
 			coverage_s = ((seq.query_end.to_f - seq.query_start.to_f)/seq.query_len.to_f) * 100
 		end
 		#if the result pass the coverage threshold
 		if coverage_s >= coverage_limit 
 			#The scrips puts another message to let us know that a good result has been found
 			puts "This sequence " + "#{query.entry_id}" + " fits the " + "#{blast_type}" + " coverage threshold."
 			#The function only returns the ID of the best hit for each blast if it passes the coverage threshold
 			return report.hits[0].definition.split("|")[0].strip
 			
 		end
  	end
 	
end


#Funtion to perform the two reciprocal blast between the two target species
#Again the input of this function is the arrays obtaining in the first funcion --> 
#these list contain for each file the file name, the sequences type that conform the file and the flat file of each proteome file 
#--> [file,f_type_,file_flat]
def do_reciprocalblasts(db1 = @file1_list, db2 = @file2_list)
	
	#information needed for both blasts 
	#the flat files for each proteome are in the third[2] position of each input array
 	flat1 = db1[2]
 	flat2 = db2[2]
 	#The names of the created databases are in the first[0] position of each input array
 	database1 = db1[0]
 	database2 = db2[0]
 	#calling 'search_blast_type' two store the types of the two blast that have to been performed in a variable
 	#for this function it is needed the same variable inputs names of the actual function 
 	blast_type1 = search_blast_type(db1, db2)
 	blast_type2 = search_blast_type(db2, db1)
 	#Create the hashes that will store the blast results
 	#As we will perform two different blast, we need two different hashes
 	@dict_hits_1 = Hash.new
 	@dict_hits_2 = Hash.new

 	#First blast
 	#We need to use all the sequences in th proteome file as queries sequences, so we have to iterate over the sequences in the flat file of the first proteome
 	#To reduce the script runung time, the first blast uses the S.pombe proteome sequences as queries and the database against which the blast is done is the one belonging to the TAIR specie
 	#Iterate over S.pombe flat file, each sequence is called query
 	flat1.each_entry do |query|
 		#for each query --> perform the blast againt Tair db, thi step is done by the 'blast' fucntion
 		#to use this script you need to change the folder path which contains your database file 
 		#this is MY PATH
 		best_hit = blast(blast_type1, query, "#{database2}")
 		#if there is not result
 		if best_hit == nil
 		#do nothing
 		 next
 		#if there is result
 		else
 			#add a new entry to the first hash, the query sequence (S.pombe) is the key, and the TAIR ID sequence (the best hit) is the value
 			@dict_hits_1[query.entry_id] = best_hit
 		end
 	end
 	#Perform the second blast: query sequences from Tair specie and database from S.pombe
 	#This time, the blast will only be performed if the Tair query sequence is a value of the first dict, to avoid doing unnecesry blasts,

 	#Iterate over Tair flat file, each sequence is called query
	flat2.each_entry do |query|
		#only if query is a value of S.pombe against Tair blast results hash
		if @dict_hits_1.values.include? query.entry_id
			#perform blast againt S.pombe database
			#to use this script you need to change the folder path which contains your database file 
			best_hit = blast(blast_type2, query, "#{database1}")
			#if there is not result
 			if best_hit == nil
 			#do nothing
 				next
 			#if there is result
 			else
 				#add a new entry to the second hash, the query sequence (Tair) is the key, and the S.Pombe ID sequence (the best hit) is the value
 				@dict_hits_2[query.entry_id] = best_hit
 			
 			end
 		#if the query sequence is not a value of the first hash
 		else
 			#do nothing
 			next
 		end
 	end
  	
  	#The function return both hashes of the two performed blasts
 	return @dict_hits_1
 	return @dict_hits_2
end

#calling the funtion to do the blasts
do_reciprocalblasts()

#And last we have to search the orthologs, so we need to compare the blast results in both hashes
#the funtion input are the two hashes created with 'do_reciprocalblasts' function
def search_orthologs(dict1 = @dict_hits_1, dict2 = @dict_hits_2)
	#in the recommended case: 
	#first dict has the following structure: dict1 --> [S.pombe]: tair
	#second dict has the following structure: dict2 --> [tair]: S.pombe
	#creating a hash to store the orthologues found between the two target species
	@orthologues = Hash.new
	#iterate over first dict, the keys are the S.pombe sequences and the values the Tair sequences
	dict1.each do |pombe, tair|
		#if the values of the dict 1 are keys in dict2; and the key searched has the sequence ID (value) that is the key in the dict1
		#in other words if the TAIR ID sequences correspond to the same S.pombe sequence in both hashes and viceversa
		 if dict2[tair] == pombe
		 	#we found an ortholog!
		 	#in this case the keys of orthologues hash are the TAIR sequences and the values the S.pombe sequences
			@orthologues[tair] = pombe
		end
	end
	#The function returns the orthologues hash 
	return @orthologues
end

#Call the function to search the orthologues
search_orthologs()

#Write the report 
#Create the report file --> csv file
File.open('orthologues_TAIR_Spombe.csv', 'w+') do |w|
	#Write an initial message to put the user in context
	w.write("Orthologue pairs between species Arabidopsis and S. pombe found by using their complete proteomes:" + "\n")
	#counter to know the final number of orthologues found 
	@contador2 = 0
	#write a header
	w.write("TAIR" + "\t" + ";" + "\t" + "S.POMBE" + "\n")
	#iterate over orthologues hash
	@orthologues.each do |tair,pombe|
		#write the othologue sequences 
		w.write("#{tair}" + "\t" + ";" + "\t" + "#{pombe}" + "\n")
		#for each pair of orthologs sum 1 to the counter
		@contador2 +=1
	end
	w.write("\n")
	#report the number of orthologs
	w.write("There are #{@contador2} orthologues between Arabidopsis and S. pombe." + "\n")
	#Stop the script execution time counter
	@ending = Process.clock_gettime(Process::CLOCK_MONOTONIC)
	#obtain the running time in seconds (time of ending - time of starting)
	@elapsed = @ending - @starting
	#obtain the running time in hours
	@elapsed_h = @elapsed.to_f / 3600
	w.write("\n")
	#write the script runtime in the report for more information
	w.write("The script runtime is " + "#{@elapsed_h} hours." + "\n")
end

#Print the number of founded orthologues and the script runtime in the terminal 
puts "There are #{@contador2} orthologues between Arabidopsis and S. pombe."
puts "The script runtime is #{@elapsed_h} hours."




