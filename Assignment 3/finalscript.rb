require 'bio'
require 'rest-client'

file = ARGV[0]

#
#Function to  access to a web page (RestClient Request) and avoid errors, this function is provided by "Bioinformatic Challenges" course.
#
# @param [String] url URL adress
# @param [Hash] headers headers, default "*/*"
# @param [String] user username for private URLs, default ""
# @param [String] pass password, default ""
#
# @return [String, false] response to URL, or false if it is something wrong 
#
def fetch(url, headers = {accept: "*/*"}, user = "", pass="")
    response = RestClient::Request.execute({
      method: :get,
      url: url.to_s,
      user: user,
      password: pass,
      headers: headers})
    return response
    
    rescue RestClient::ExceptionWithResponse => e
      $stderr.puts e.inspect
      response = false
      return response  # now we are returning 'False', and we will check that with an \"if\" statement in our main code
    rescue RestClient::Exception => e
      $stderr.puts e.inspect
      response = false
      return response  # now we are returning 'False', and we will check that with an \"if\" statement in our main code
    rescue Exception => e
      $stderr.puts e.inspect
      response = false
      return response  # now we are returning 'False', and we will check that with an \"if\" statement in our main code
end

#
#Function which reads the "ArabidopsisSubNetwork_GeneList.txt", which contains more than 100 genes in in AGI locus code format (one in each column),
#creates an array (genes) containin all the file genes in lowercase letters (to make easier the id search).
#Then the function search the EMBL entry for each gene of the array (genes) and creates a hash to store the gene id as key and its embl entry as value.
#In case that there were no response for the gene URL the hash would record the gene id as key and the message "There is no respone for the gene (gene id)" as value.
#
# @param file [String] path to genes_file of target genes, default "command line (possition [0])"
#
# @return [Array<Symbol>] an Array (genes) of gene ids of the input file
# @return [Hash<Symbol => Bio::EMBL, String>] a Hash (@gene_page_hash) containing the genes id and its embl entries as a Bio::EMBL object or a message if there is no response from the URL
#

def search_gene_web(file) 
  genes = Array.new #create a new array
  File.open(file).each_line do |gene| #open file and go through the file line by line 
      gene = gene.strip #to no create lists and delete line breaks
      gene = gene.downcase #(downcase) to match with the results of the web page
      genes << gene #add each id to the genes array
  end 

  @gene_page_hash = Hash.new #create a new hash
  genes.each do |id| #scroll through the list of genes
    #search each id in ensembl with the URL
    response = fetch("http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=ensemblgenomesgene&format=embl&id=#{id}") # this is the URL of ensembl
    
    #If there is response from the URL record the entry of each gene as a Bio::EMBL object
    if response
      embl_gen = Bio::EMBL.new(response.body) #create the Bio::EMBL object
      @gene_page_hash[id] = embl_gen #record the gene id and the Bio::EMBL object in the hash

    #if there is no response from the URL, record the id and a message in the hash 
    else
      @gene_page_hash[id] = ["There is no web response for the gene #{id}"]
    end
  end
  return genes 
  return @gene_page_hash
end

#Call the first function
search_gene_web(file)
#
#Funtion to search the gene exons of a given hash (@gene_page_hash) [gene => Bio::EMBL, String], search the sequence "CTTCTT" in all exons and record the start position of all the coincidences in two different arrays, 
#one for the positive strand and the other for the negative one (for each gene).
#The positions, are first searched for the exon, and then located in the gene absolute sequence. 
#The function returns a hash [gene => [[Bio::EMBL], [position_forward],[position_reverse]]].
#For the forward strand the target secuence is cttctt and for the reverse strand the target sequence is its complementary reverse (aagaag).
#It is compulsory to sum 1 to the positions because ruby return the coincidence initializing in the position 0 not 1, and it is required the real genomic position. 
#The function also creates an array for those genes that do not have any sequence_cttctt coincidence in their exons.
#
# @param hash_wb [Hash<Symbol => Bio::EMBL, String>] a Hash, default @gene_page_hash
#
# @return [Hash<Symbol => Array<Array<Bio::EMBL, integrer, integrer>>] hash (@hash_pos) whose keys are again the gene id, and the value is an array of arrays: Bio::EMBL object, array for forward positions and array for reverse positions
# @return [Array<Symbol>] array (@lista_no) of gene ids without cttctt sequence in its exons
#

def search_repeats(hash_wb = @gene_page_hash)
  @hash_pos = Hash.new #Create a new hash to obtain forward and reverse repeat (cttctt) start positions for each gene
  #scroll the hash containing the genes and their Bio::EMBL objects
  hash_wb.each do |gene_id, embl_entry|
    #create an array for the forward positions of each gene
    list_posf = []
        
    ##create an array for the reverse positions of each gene
    list_posr = []

    #search gene sequence using BioRuby
    gen_seq = embl_entry.seq
  
    #search in the object features and filter by exon category because we are looking for repeats in exons (BioRuby)
    embl_entry.features.each do |feature| 
      if feature.feature == 'exon'
        #search the location coords. of the exon 
        feature.locations.each do |location|
          #store this location in a varible becuase it is easy to organize myself this way
          coor = location
          #obtain the exon sequence 
          exon_seq = gen_seq[coor.from..coor.to]
          
          #if there is no sequence for the exon, do not enter in the loop
          if exon_seq.nil? == false

          #if we are in the positive strand search 'cttctt' sequence
            if location.strand == 1
              #In order to search all the start position of the repetition ('cttctt') coincidences, inlcuding those which are overlapping, I searched this function:
              #https://stackoverflow.com/questions/43329481/find-all-indices-of-a-substring-within-a-string
              #'cttctt' is the Regexp searched 
              positions = exon_seq.enum_for(:scan, /(?=cttctt)/).map { Regexp.last_match.offset(0).first}
              #if there is any 'cttctt' coincidence in the exon, skip 
              if positions.empty? == false
                #transform matching positions in the exon to the position in the gene and add 1 (the first position should be counted as 1, not 0)
                positions_f = positions.map {|pos| [((pos +  coor.from) + 1)]}
                #only include the founded positions to the forward_positions_array if they are not included yet, to avoid possible repetitions
                next if list_posf.include?(positions_f)
                list_posf << positions_f

              end
            
            #if we are in the negative strand search 'aagaag' sequence, in other words, 'cttctt' complementary reverse
            elsif location.strand == -1
              #In order to search all the start position of the repetition ('aagaag') coincidences, inlcuding those which are overlapping, I searched this function:
              #https://stackoverflow.com/questions/43329481/find-all-indices-of-a-substring-within-a-string
              #'aagaag' is the Regexp searched 
              positions2 = exon_seq.enum_for(:scan, /(?=aagaag)/).map { Regexp.last_match.offset(0).first}
              #if there is any 'aagaag' coincidence in the exon, skip 
              if positions2.empty? == false
                #transform matching positions in the exon to the position in the gene and add 1 (the first position should be counted as 1, not 0)
                positions_r = positions2.map {|pos| [((pos +  coor.from) + 1)]}
                #only include the founded positions to the reverse_positions_array if they are not included yet, to avoid possible repetitions 
                next if list_posr.include?(positions_r)
                list_posr << positions_r
              end
            end
          end
        end
      end
    end
    #when both positive and negative string matches have already been searched for a gene and you have an array of arrays for each of the strings,
    #these arrays of arrays are transformed in a simple array of start matching positions, in order to descard the repeated positons that appear when we transform exon positions to gen positions
    if !list_posr.empty? == true
      #obtain a simple array
      list_posr = list_posr.flatten
      #delete repeated postions
      list_posr = list_posr.uniq
    end

    if !list_posf.empty? == true
      list_posf = list_posf.flatten
      list_posf = list_posf.uniq
    end
    #create a hash entry for each gene, the value is an array of arrays: position 0 is the Bio::EMBL object, position 1 is the array for the start matching postions in the positive strand 
    #and position 2 is the the array for the start matching postions in the negative strand
    @hash_pos[gene_id] = [embl_entry, list_posf, list_posr]
    
    #create an array for those genes that do NOT have exons with the 'CTTCTT' repeat
    @lista_no = []
    #scroll the hash containing the positions where the target sequence starts for both strands and for all genes
    @hash_pos.each do |key, value|
      #if both, the array of forward and reverse positions are empty, this means that there is no match 
      if value[1].empty? == true and value[2].empty? == true
        #include the gene in the array 
        @lista_no << key
      end
    end
  end
  return @hash_pos
  return @lista_no 
end

search_repeats()

#
#Funtion which receives a hash (as input) with the gene_ids as keys and lists of lists as values: Bio::Embl object[0] and their start positions (forward [1] and reverse[2]) for the 'cttctt' sequence,
#and create a Bio::sequence object from the Bio::EMBL object for each gene, thus, it is possible to add new features to the object.
#Then, the function adds new 'myrepeat' features for each cttctt sequence found in each gene exon.  
#The new feature contains the start and end position of the cttctt sequence in the gene (position_), repeat_motif ('cttctt' sequence), gene id (gene), strand, and an ID.
#This function returns a hash [gene => Bio::seq] with the new "myrepeat" features for each gene. 
#The function "knows" if the target sequence is in the positive or negative strand depending on the position of the hash values that receives: positive, if list[1] and negative, if list[2].
#
# @param info [Hash<Symbol => Array<Array<Bio::EMBL, integrer, integrer>>] a Hash (@hash_pos) whose keys are again the gene id, and the value is an array of arrays: Bio::EMBL object, array for forward positions and array for reverse positions, default @hash_pos
#
# @return [Hash<Symbol => Bio::Sequence>] a Hash (@total_hash) whose keys are the gene ids and the values its Bio::Sequence objects with the new 'myrepeat' feature for those who have the 'cttctt' sequence in their exons
# 

def new_feature1(info = @hash_pos)
  #create a new hash to record the gene id and its Bio::Sequence object with all its features 
  @total_hash = Hash.new
  #scroll the  hash containig the gene id and an array of arrays: position 0 is the Bio::EMBL object, position 1 is the array for the start matching postions in the positive strand 
  #and position 2 is the the array for the start matching postions in the negative strand, for each gene  
  info.each do |gene, lista| 

    #create a Bio::Sequence object (BioRuby) from the Bio:: EMBL object to be able to add new features
    bioseq_obj = lista[0].to_biosequence
    #positions on positive strand are in the second array position
    positive_rep = lista[1]
    #positions on negative strand are in the third array position
    negative_rep = lista[2]

    #go through all the positions of the positive strand
    positive_rep.each do |pos2|
      #search the index of the position to obtain an unique ID for each matching position and sum 1 to not start in 0
      indice = (positive_rep.index(pos2).to_i) + 1
      #store the start position
      start_position = pos2
      #obtain the end position by sum the length - 1 of the sequence, because the end position is the position of the last T('CTTCTT' = 6 -1 = 5)
      end_position = pos2.to_i + 5 
      #obtain the gene id in upcase letters 
      gene_s = gene.upcase
      #obtain the position in the format that appears in the EMBL web
      position_ = "#{start_position}..#{end_position}"
      #if were are in the array of positive strand, the strand symbol is '+'      
      strand = "+"
      #create a new Bio::Sequence object feature called 'myrepeat' and the position (start and end) of the repetition in the gene, as it is done in the Jupyter Notebooks of the 'Challenges Course'
      seq_feature = Bio::Feature.new("myrepeat", position_)
      #append the repeat motif ('CTTCTT') to the feature 
      seq_feature.append(Bio::Feature::Qualifier.new("repeat_motif", "CTTCTT")) 
      #append the gene id to the feature 
      seq_feature.append(Bio::Feature::Qualifier.new("gene",gene_s))
      #append the strand value ('+') to the feature 
      seq_feature.append(Bio::Feature::Qualifier.new("strand", strand))
      #append the unique ID value for each matching to the feature 
      seq_feature.append(Bio::Feature::Qualifier.new("ID", "ID =" + "found_exon_repetition=CTTCTT_in" + "#{gene_s}" + "." + "#{indice}" + "/." + "#{strand}." + ";"))
      #add the new feature to the Bio::Sequence object
      bioseq_obj.features << seq_feature
    end


    #go through all the positions of the negative strand
    negative_rep.each do |pos2|
      #search the index of the position to obtain an unique ID for each matching position and sum 1 to not start in 0
      indice = (negative_rep.index(pos2).to_i) + 1
      #store the start position
      start_position = pos2
      #obtain the end position by sum the length - 1 of the sequence, because the end position is the position of the last T('AAGAAG' = 6 -1 = 5)
      end_position = pos2.to_i + 5    
      #obtain the gene id in upcase letters     
      gene_s = gene.upcase
      #obtain the position in the format that appears in the EMBL web and add "complement" string as we are in the minus strand
      position = "complement(#{start_position}..#{end_position})"
      #if were are in the array of negative strand, the strand symbol is '-' 
      strand = "-"
      #create a new Bio::Sequence object feature called 'myrepeat' and the position (start and end) of the repetition in the gene, as it is done in the Jupyter Notebooks of the 'Challenges Course'
      seq_feature = Bio::Feature.new("myrepeat", position)
      #append the repeat motif ('CTTCTT') to the feature
      seq_feature.append(Bio::Feature::Qualifier.new("repeat_motif", "CTTCTT"))
      #append the gene id to the feature
      seq_feature.append(Bio::Feature::Qualifier.new("gene", gene_s))
      #append the strand value ('-') to the feature
      seq_feature.append(Bio::Feature::Qualifier.new("strand", strand))
      #append the unique ID value for each matching to the feature 
      seq_feature.append(Bio::Feature::Qualifier.new("ID", "ID =" + "found_exon_repetition=CTTCTT_in" + "#{gene_s}" + "." + "#{indice}" + "/." + "#{strand}." + ";"))
      #add the new feature to the Bio::Sequence object
      bioseq_obj.features << seq_feature
      
    end
    
    #record all the genes id and their Bio::Sequence object with the new features 
    @total_hash[gene] = bioseq_obj
  end
  return @total_hash
end

new_feature1()


#
# Funtion to create a GFF3-formatted file with all of the new 'myrepeat' features for each gene.
# The GFF3 fields are the following; sequence ID, source, feature type, feature start, feature end, score, strand, phase, atributes. 
#
# @param info_hash [Hash<Symbol , Bio::Sequence>]  a Hash of Bio::Sequence objects (containing 'myrepeat' features), default @total_hash
#
# @return [String] message for knowing that the function was succesfully called
# 

def write_gff3(info_hash = @total_hash)
  #create a new gff3 format file to be written 
  gff3 = File.new("gene_positions.gff3", "w")
  #write the first line of the file 
  gff3.write("##gff-version 3\n")
  #The gff3 file will have 9 different columns;
  # A source: describes the algorithm or the procedure that generated this feature, I decided to write the page from where i obtained de informtion (EMBL)
  source = "EMBL"
  # type: describes what the feature is (mRNA, domain, exon, etc.), it is a repeated sequence found in a gene exon 
  type = "direct_repeat_in_EXON"
  # score: typically E-values for sequence similarity and P-values for predictions, its values is '.'
  score = "."
  #phase: indicates where the feature begins with reference to the reading frame, its values is '.'
  phase = "."
  #scroll the hash containing the genes id and their Bio::Sequence object (with the new features)
  info_hash.each do |g, object|
    #search within the object's features 
    object.features.each do |feature|
      featurename = feature.feature
      #and filter by the feature 'myrepeat'
      next unless featurename == "myrepeat"
      #search (BioRuby) the location of the feature (position of the sequence matched in the gene)
      location_final = feature.locations[0]
      #record star postion (BioRuby)
      start_ = location_final.from
      #record the end position (BioRuby)
      end_ = location_final.to
      #call the associated categories of the 'myrepeat' feature (BioRuby)
      #call gene id and record it as the sequence ID for the gff3 file
      seqid = (feature.assoc["gene"])
      #call ans record strand for the unique ID and strand fields of the gff3 file 
      strand = (feature.assoc["strand"]).to_s
      #call and record the ID for the gff3 file 
      attributes = (feature.assoc["ID"]).to_s
      #write the gff3 lines: one for each matching of 'CTTCTT'repeat in the exon sequence (all fields must be separated by '\t')
      gff3.write("#{seqid}" + "\t" + "#{source}" + "\t" + "#{type}" + "\t" + "#{start_}" + "\t" + "#{end_}" + "\t" + "#{score}" + "\t" + "#{strand}" + "\t" + "#{phase}" + "\t" + "#{attributes}" + "\n")
    end
  end
  return "The gene_positions.gff3 file was created "

end

#call the function for writing the gff3 file 
write_gff3()

#
# Funtion to write a report showing which (if any) genes on the target list do NOT have exons with the 'CTTCTT' repeat. 
#
# @param genes_no [Array<Symbol>] an Array (@lista_no) of gene ids without 'cttctt' sequence in its exons
#
# @return [String] message for knowing that the function was succesfully called
#

def write_nogenes (genes_no = @lista_no)
  #create a new txt format file to be written 
  txt = File.new("norepeats.txt", "w")
  txt.write("\n")
  txt.write("\n")
  #Write "REPORT"
  txt.write("\t" + "\t" + "\t" + "\t" + "\t" + "REPORT" + "\n")
  txt.write("\n")
  txt.write("\n")
  #write a description
  txt.write("This file is for reporting those genes of the list that do NOT have exons with the CTTCTT repeat". + "\n")
  txt.write("\n")
  #record the number of genes that do NOT have exons with the 'CTTCTT' repeat
  num_gen = genes_no.length()
  #write the number of genes that do NOT have exons with the 'CTTCTT' repeat
  txt.write("There are " + "#{num_gen}" + " genes without exons with the CTTCTT repeat." + "\n")
  txt.write("\n")
  #scroll the array containing all the genes that do NOT have exons with the 'CTTCTT' repeat
  genes_no.each do |gene|
    #write each gene of the array
    txt.write("\t" + "\t" + "- " + "#{gene}" + "\n")
    txt.write("\n")
  end
  return "The norepeats.txt file was created "
end

#call the function which writes the report file 
write_nogenes()

#
#Funtion which receives a hash (as input) with the gene_ids as keys and lists of lists as values: gen id Bio::Embl object[0] and its start positions (forward [1] and reverse[2]) for the 'cttctt' sequence,
#and create a Bio::sequence object from the Bio::EMBL object for each gene, thus, it is possible to add new features to the object.
#Then, the function adds new 'myrepeat_chr' features for each cttctt sequence found in each gene exon, in this case the position is not for the gene, is for full chromosome coordinates used by EnsEMBL.  
#This function returns a hash [gene => Bio::seq] with the new "myrepeat_chr" features for each gene. 
#The new feature contains the start and end position of the cttctt sequence in the chromosome (position), chromosome number (chr), repeat_motif ('cttctt' sequence), gene id (gene), strand, and an ID.
#The function "knows" if the target sequence is in the positive or negative strand depending on the position of the hash values that receives: positive, if list[1] and negative, if list[2].
#
# @param info2 hash [Hash<Symbol => <Array<Array<Bio::EMBL, integrer, integrer>>>]  a Hash (@hash_pos) whose keys are again the gene id, and the value is an array of arrays: Bio::EMBL object, array for forward positions and array for reverse positions, default @hash_pos
#
# @return [Hash<Symbol => Bio::Sequence>] hash (@total_hash2) whose keys are the gene ids and the values its Bio::Sequence objects with the new 'myrepeat_chr' feature for those who have the 'cttctt' sequence in their exons
# 

def chr(info2 = @hash_pos)
  #create a new hash to record the gene id and its Bio::Sequence object with all its features 
  @total_hash2 = Hash.new
  #scroll the  hash containig the gene id and an array of arrays: position 0 is the Bio::EMBL object, position 1 is the array for the start matching postions in the positive strand 
  #and position 2 is the the array for the start matching postions in the negative strand, for each gene  
  info2.each do |gen_i, list_rep| 
    #create a Bio::Sequence object (BioRuby) from the Bio:: EMBL object to be able to add new features, in this case the searched position is not for the gene, is for the whole chromosome 
    entry = list_rep[0].to_biosequence
    #searcg the chromosome number in the Bio::Sequence object(BioRuby)
    chr_num = entry.primary_accession.split(":")[2]
    #searcg the chromosome start position in the Bio::Sequence object(BioRuby)
    chr_to = entry.primary_accession.split(":")[3]
    #positions on positive strand are in the second array position
    positive_rep = list_rep[1]
    #positions on negative strand are in the third array position
    negative_rep = list_rep[2]

    #go through all the positions of the positive strand
    positive_rep.each do |position|
      #search the index of the position to obtain an unique ID for each matching position and sum 1 to not start in 0
      indice = positive_rep.index(position)
      #store the start position in the gen and transfom it into the start postion in the chromosome 
      start_pos = (chr_to.to_i) +  (position).to_i - 1 
      #obtain the end position by sum the length - 1 of the sequence, because the end position is the position of the last T('CTTCTT' = 6 -1 = 5)
      end_pos = start_pos.to_i + 5
      #obtain the gene id in upcase letters 
      gene_s = gen_i.upcase
      #obtain the position in the format that appears in the EMBL web 
      position = "#{start_pos}..#{end_pos}"
      #if were are in the array of positive strand, the strand symbol is '+' 
      strand = "+"
      #create a new Bio::Sequence object feature called 'myrepeat_chr' and the position (start and end) of the repetition in the chromosome, as it is done in the Jupyter Notebooks of the 'Challenges Course'
      seq_feature2 = Bio::Feature.new("myrepeat_chr", position)
      #append the number of the chromosome (chr_num)
      seq_feature2.append(Bio::Feature::Qualifier.new("chr", "chr" + "#{chr_num}"))
      #append to the feature the repeat motif ('CTTCTT') to the feature
      seq_feature2.append(Bio::Feature::Qualifier.new("repeat_motif", 'CTTCTT'))
      #append the gene id to the feature
      seq_feature2.append(Bio::Feature::Qualifier.new("gene", gene_s))
      #append the strand to the feature ('+')
      seq_feature2.append(Bio::Feature::Qualifier.new("strand", strand))
      #append the unique ID value for each matching to the feature 
      seq_feature2.append(Bio::Feature::Qualifier.new("ID", "ID =" + "#{gene_s}" + "." + "#{chr_num}" + "." + "#{(indice + 1)}" + "found_exon_repetition=CTTCTT" + "/." + "#{strand}" + ";"))
      #add the new feature to the Bio::Sequence object
      entry.features << seq_feature2
    end
  

    #go through all the positions of the negative strand
    negative_rep.each do |position|
      #search the index of the position to obtain an unique ID for each matching position and sum 1 to not start in 0
      indice = negative_rep.index(position)
      #store the start position in the gen and transfom it into the start postion in the chromosome 
      start_pos = (chr_to.to_i) + (position).to_i - 1 
      #obtain the end position by sum the length - 1 of the sequence, because the end position is the position of the last T('AAGAAG' = 6 -1 = 5)
      end_pos = start_pos.to_i + 5
      #obtain the gene id in upcase letters 
      gene_s = gen_i.upcase
      #obtain the position in the format that appears in the EMBL web and add "complement" string as we are in the minus strand
      position = "complement(#{start_pos}..#{end_pos})"
      #if were are in the array of negative strand, the strand symbol is '-' 
      strand = "-"
      #create a new Bio::Sequence object feature called 'myrepeat_chr' and the position (start and end) of the repetition in the chromosome, as it is done in the Jupyter Notebooks of the 'Challenges Course'
      seq_feature2 = Bio::Feature.new("myrepeat_chr", position)
      #append the number of the chromosome (chr_num)
      seq_feature2.append(Bio::Feature::Qualifier.new("chr", "chr" + "#{chr_num}"))
      #append to the feature the repeat motif ('CTTCTT') to the feature
      seq_feature2.append(Bio::Feature::Qualifier.new("repeat_motif", 'CTTCTT'))
      #append the gene id to the feature
      seq_feature2.append(Bio::Feature::Qualifier.new("gene", gene_s))
      #append the strand to the feature ('-')
      seq_feature2.append(Bio::Feature::Qualifier.new("strand", strand))
      #append the unique ID value for each matching to the feature
      seq_feature2.append(Bio::Feature::Qualifier.new("ID", "ID =" + "#{gene_s}" + "." + "#{chr_num}" + "." + "#{(indice + 1)}" + "found_exon_repetition=CTTCTT" + "/." + "#{strand}" + ";"))
      #add the new feature to the Bio::Sequence object
      entry.features << seq_feature2
    end
    #record all the genes id and their Bio::Sequence object with the new chr features
    @total_hash2[gen_i] = entry
  end
  return @total_hash2
end

chr()


#
# Funtion to create a GFF3-formatted file with all of the new 'myrepeat_chr' features for each gene.
# The GFF3 fields are the following; sequence ID (chromosome), source, feature type, feature start, feature end, score, strand, phase, atributes. 
#
# @param info_chr [Hash<Symbol, Bio::Sequence>] a Hash of Bio::Sequence objects, containing 'myrepeat_chr' features, default @total_hash2
#
# @return [String] message for knowing that the function was succesfully called
#

def write_gff3_chr(info_chr = @total_hash2)
  #create a new gff3 format file to be written
  #fields are the same that for the gene position gff3 file, but in this case the sequence ID is the chomosome number and the positions are chromosomal positions
  gff3_chr = File.new("chr.gff3", "w")
  gff3_chr.write("##gff-version 3\n")
  source = "EMBL"
  type = "direct_repeat_in_EXON"
  score = "."
  phase = "."
  #scroll the hash containing the genes id and their Bio::Sequence object (with the new chr features)
  info_chr.each do |g, object|
    
    object.features.each do |feature|
      featuretype = feature.feature
      #and filter by the feature 'myrepeat_chr'
      next unless featuretype == "myrepeat_chr"
      location2 = feature.locations[0]
      start = location2.from
      end_ = location2.to
      strand = (feature.assoc["strand"]).to_s
      attributes = (feature.assoc["ID"]).to_s
      seqid = (feature.assoc["chr"]).to_s
      gff3_chr.write("#{seqid}" + "\t" + "#{source}" + "\t" + "#{type}" + "\t" + "#{start}" + "\t" + "#{end_}" + "\t" + "#{score}" + "\t" + "#{strand}" + "\t" + "#{phase}" + "\t" + "#{attributes}" + "\n")
    end
  end
  return "The chr.gff3 file was created "
end

#call the function for writing the gff3 file
write_gff3_chr()

  

