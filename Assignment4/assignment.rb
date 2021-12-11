require 'bio'

filearab = ARGV[0]
filespombe = ARGV[1]

def important_info (file1 = filearab, file2 = filespombe)
	@pombe_list = Array.new
	@arab_list = Array.new
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
	
		file_out = (file.to_str).chomp(".fa") + "_out.fa"
		system("makeblastdb -in #{file} -dbtype '#{f_type_}' -out #{file_out}")
		contador +=1
		if contador == 1
			@arab_list = [file,file_out,f_type_,file_flat]

		else
			@pombe_list = [file,file_out,f_type_,file_flat]
	
		end
	end
	return @arab_list
	return @pombe_list

end

important_info(filearab,filespombe)

### Ya tenemos las bases de datos de los respectivos archivos
#Creamos una factoria dependiendo de en qué dirección se hace el blast 
#Si es de arab contra pombe -->
def search_blast_type(data1 = @arab_list, data2 = @pombe_list)
	type1 = data1[2]
	type2 = data2[2]

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

#hacer un blast dependiendo de lo que se quiera le daremos un valor a db1 u otro 

#queremos llamar a esta función dos veces para hacer un blast recíproco 
def blast_db1_db2(db1 = @arab_list, db2= @pombe_list)
	best_hit_dict = Hash.new
	blast_type = search_blast_type(db1, db2)
	db_blast = db2[1]
	factory = Bio::Blast.local("#{blast_type}", "#{db_blast}","-F 'm S' -e 1e-6") 
	
	#como las secuencias que estamos blasteando son de la db1 pues buscamos en el flat de la db1 (en la posición 3 de la lista )
	flat_db1 = db1[3]
	flat_db1.each_entry do |query|
		report = factory.query(query)
        if !report.hits.empty? 
        	#esta se saca del report del blast
       	 	seq_procedente_db2 = report.hits.first.definition.split("|")[0].split
        	best_hits_dict[query.entry_id] = seq_procedente_db2 
    end
    return best_hits_dict
 end

#Ahora para buscar los ortólogos se comparan los resutados del blast en ambas direcciones 

def get_best_reciprocal_hits(db1_hash = @arab_list, db2_hash = @pombe_list) 

    # Get the best hits from blasting every sequence of db1 against db2
    best_hits_tair_against_pombe = blast_db1_db2(@arab_list, @pombe_list) 
    # Get the best hits from blasting every sequence of db2 against db1, including the previous results to avoid doing unnecessary blasts
    best_hits_pombe_against_tair = blast_db_against_db(@pombe_list, @arab_list)

    # búsqueda de ortólogos 
    @orthologs = {}
    #búsqueda en las llaves del diccionario de spombe contra tair, las llaves son genes de pombe y las seq de tair
    best_hits_pombe_against_tair.each_key do |key|
        # la secuencia de tair es el value de las llaves 
        seq1 = best_hits_pombe_against_tair[key]
        # si la secuencia de pombe coincide con el valor de la llave de tair, o sea el pombe para el diccionario de tair vs pombe 
        if key == best_hits_q1_against_db2[seq1] 
        	#tenemos ortólogo
        	@orthologs[key] = seq1 # storing the results on the hash
    end
    return reciprocal_hits 
end

get_best_reciprocal_hits()

File.open('orthologues_TAIR_Spombe.txt', 'w+') do |w|
	w.write("Orthologue pairs between species Arabidopsis and S. pombe found by using their complete proteomes:")
	contador2 = 0
	@orthologs.each do |seq1,seq2|
		w.write("TAIR SEQ" + "\t" + ":" + "\t" + "S.POMBE")
		w.write("#{seq1}" + "\t" + ":" + "\t" + "#{seq2}")
		contador2 +=1
	end
	w.write("There are #{contador2} orthologues between Arabidopsis and S. pombe")
end



