#Link classes files to main script

require "./seedstock.rb"
require "./crossclass.rb"
require "./geneclass.rb"

#declare input files
seed_stock_file = ARGV[0]
cross_data_file = ARGV[1]
gene_information_file = ARGV[2]

#First task
s = Seedstock.new()

#Define "contador" variable for storaging header line of seed_stock_file
contador1 = 1

#Create an array to store the different s objects
array_stock = []

#Read seed_stock_file and storage the elements of each line as object properties
File.readlines(seed_stock_file).each do |line| 
  #Header is a constant that storage the header line of the file
  if contador1 == 1
    @header = line
    contador1 += 1
    
  else
    line = line.split("\t")
    seedrem = Integer(line[4]) - 7
    s = Seedstock.new(
      :seed_stock => line[0], 
      :mutant_gene_id => line[1], 
      :last_planted => line[2], 
      :storage => line[3], 
      :grams_remaining => Integer(line[4])
      )
    array_stock = array_stock.push(s)
  end
end


#Create and opening a new .tsv file 
File.open("seed_stock_data_2.tsv", "w") {|a|
  a.write @header
  #Go down the array_stock for write all objects
  array_stock.each do |elementt|
    i = array_stock.find_index(elementt)
    #Call methods
    #In this case the number given to actualize_rem_seed_and_date method is seven because we want to plant 7 grams
    array_stock[i].actualize_rem_seed_and_date(7)  
    #This method writes all properties in the new file
    a.write array_stock[i].write_file
  end
  }

###########################################################
c = Crossdata.new()

#Read cross_data_file and storage the elements of each line as object properties
contador2 = 1
array = []
File.readlines(cross_data_file).drop(1).each do |line| 
  line = line.split("\t")
  c = Crossdata.new(
    :parent1 => line[0], 
    :parent2 => line[1], 
    :f2_wild => line[2], 
    :f2_p1 => line[3], 
    :f2_p2 => line[4],
    :f2_p1p2 => line[5],
    )
    
  array = array.push(c)

end

#Calculate chi for all the c objects and put Final report messages 
array.each do |element|
  i = array.find_index(element)
  array[i].chi_cuadrado

  
end

###########################################################

#Read gene_information_file and storage the elements of each line as object properties

array_gene = []
File.readlines(gene_information_file).drop(1).each do |line| 
  line = line.split("\t")
  g = Geneinformation.new(
    :gene_id => line[0], 
    :gene_name => line[1], 
    :mutant_phenotype => line[2], 
    )
    
  array_gene = array_gene.push(g)

end
