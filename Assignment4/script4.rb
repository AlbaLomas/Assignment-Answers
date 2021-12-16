require 'bio'

@starting = Process.clock_gettime(Process::CLOCK_MONOTONIC)

file1 = ARGV[0]
file2 = ARGV[1]

def important_info (file1 = file1, file2 = file2)
	@file1_list = Array.new
	@file2_list = Array.new

	files = [file1, file2]
	contador = 0
	files.each do |file|

		file_flat = Bio::FlatFile.auto(file)
		file_type = file_flat.next_entry.to_biosequence.guess
		if file_type == Bio::Sequence::NA
			f_type_ = 'nucl'
		elsif file_type == Bio::Sequence::AA
			f_type_ = 'prot' 
		end
	
		system("makeblastdb -in #{file} -dbtype '#{f_type_}'")
	
		contador +=1
		if contador == 1
			@file1_list = [file,f_type_,file_flat]

		else
			@file2_list = [file,f_type_,file_flat]
	
		end
	end
	return @file1_list
	return @file2_list

end

important_info(file1,file2)

#search blast type

def search_blast_type(data1 = @file1_list, data2 = @file2_list)
	type1 = data1[1]
	type2 = data2[1]

	if type1 == "nucl" && type2 == "prot"
		blastype = "blastx"
	elsif type1 == "prot" && type2 == "nucl"
		blastype = "tblastn"
	elsif type1 == "prot" && type2 == "prot"
		blastype = "blastp"
	elsif type1 == "nucl" && type2 == "nucl"
		blastype = "blastn"
	end
	return blastype
end

#doing blst
#para correr un blast se necesita una secuencia query, una base de datos y el tipo de blast
def blast(blast_type = nil, query = nil , db = nil)
	coverage_limit = 50
	factory = Bio::Blast.local("#{blast_type}", "#{db}","-F T -s F  -e 1e-6") 
	#como las secuencias que estamos blasteando son de la db1 pues buscamos en el flat de la db1 (en la posición 3 de la lista )
	report = factory.query(query)
	#esta se saca del report del blast
	if report.hits[0] == nil 
		puts "Report is empty for the query secuence " + "#{query.entry_id}" + " when running a " + "#{blast_type}" + " against " + "#{db.split("/")[5]}" + " database."
	else
		seq = report.hits[0]
		coverage = ((seq.query_end.to_f - seq.query_start.to_f)/seq.query_len.to_f) * 100
 		if coverage >= coverage_limit
 			return report.hits[0].definition.split("|")[0].strip
 		end
  	end
 	
end


def do_reciprocalblasts(db1 = @file1_list, db2 = @file2_list)
	
	#information need for both blasts 
 	flat1 = db1[2]
 	flat2 = db2[2]
 	database1 = db1[0]
 	database2 = db2[0]
 	blast_type1 = search_blast_type(db1, db2)
 	blast_type2 = search_blast_type(db2, db1)
 	@dict_hits_1 = Hash.new
 	@dict_hits_2 = Hash.new
 	##primer blast 
 	flat1.each_entry do |query|
 		best_hit = blast(blast_type1, query, "/Users/albalomasredondo/Desktop/pruebas/"+ "#{database2}")
 		if best_hit == nil
 		 next
 		else
 			@dict_hits_1[query.entry_id] = best_hit
 		end
 	end
 	##segundo blast pero filtrando por la secuencias de tair que se hayan encontrado en el primer blast
	flat2.each_entry do |query|
		if @dict_hits_1.values.include? query.entry_id
			best_hit = blast(blast_type2, query,"/Users/albalomasredondo/Desktop/pruebas/"+ "#{database1}")
 			if best_hit == nil
 		 		next
 			else
 				@dict_hits_2[query.entry_id] = best_hit
 			
 			end
 		else
 			next
 		end
 	end
  
 	return @dict_hits_1
 	return @dict_hits_2
end

do_reciprocalblasts()

def search_orthologs(dict1 = @dict_hits_1, dict2 = @dict_hits_2)
	#dict1 --> [S.pombe]: tair
	#dict2 --> [tair]: S.pombe
	@orthologues = Hash.new
	dict1.each do |pombe, tair|
		 if dict2[tair] == pombe
			@orthologues[tair] = pombe
		end
	end
	return @orthologues
end

search_orthologs()


File.open('orthologues_TAIR_Spombe.csv', 'w+') do |w|
	w.write("Orthologue pairs between species Arabidopsis and S. pombe found by using their complete proteomes:" + "\n")
	@contador2 = 0
	w.write("TAIR" + "\t" + ";" + "\t" + "S.POMBE" + "\n")
	@orthologues.each do |tair,pombe|
		w.write("#{tair}" + "\t" + ";" + "\t" + "#{pombe}" + "\n")
		@contador2 +=1
	end
	w.write("\n")
	w.write("There are #{@contador2} orthologues between Arabidopsis and S. pombe." + "\n")
	@ending = Process.clock_gettime(Process::CLOCK_MONOTONIC)
	@elapsed = @ending - @starting
	@elapsed_h = @elapsed.to_f / 3600
	w.write("\n")
	w.write("The script runtime is " + "#{@elapsed_h} hours." + "\n")
end

puts "There are #{@contador2} orthologues between Arabidopsis and S. pombe."
puts "The script runtime is #{@elapsed_h} hours."




