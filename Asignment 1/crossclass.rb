#class cross_data

class Crossdata 
  attr_accessor :parent1
  attr_accessor :parent2
  attr_accessor :f2_wild
  attr_accessor :f2_p1
  attr_accessor :f2_p2
  attr_accessor :f2_p1p2 
  attr_accessor :chi_2
  attr_accessor :linked
  
  
  def initialize (params = {}) # get a name from the "new" call, or set a default
    @parent1= params.fetch(:parent1, "0000")
    @parent2= params.fetch(:parent2, "00000")
    @gene_name1 = params.fetch(:gene_name1, "0000")
    @gene_name2 = params.fetch(:gene_name2, "0000")
    @f2_wild = params.fetch(:f2_wild, "000")
    @f2_p1 = params.fetch(:f2_p1, "00")
    @f2_p2 = params.fetch(:f2_p2, "00")
    @f2_p1p2 = params.fetch(:f2_p1p2, "00")
  	@chi_2 = params.fetch(:chi_2, "000.00000000000")
    @linked = params.fetch(:linked, "No linked")
  end

  def test_array()
    puts "#{@parent1}\t#{@parent2}\t#{@f2_wild}\t#{@f2_p1}\t#{@f2_p2}\t#{@f2_p1p2}\t#{@linked}"
  end

  #Calculate chi cuadrado and create hashes for connecting classes to get the final report
  def chi_cuadrado()
    linked_array = []
    #hash_name --> dict key (gene id), values (gene name)
    hash_name = {}
    name_array = []

    File.readlines("./gene_information.tsv").drop(1).each do |line|
      line = line.split("\t")
      hash_name[line[0]] = line[1]
      name_array << line[1]
    end
    #hash_id --> dict key (seed stock), values (gene id)
    hash_id = {}
    stock_array = []

    File.readlines("./seed_stock_data.tsv").drop(1).each do |line|
      line = line.split("\t")
      hash_id[line[0]] = line[1]
      stock_array << line[0]
    end

    #hash_stock_name --> dict key (seed stock), values (gene name)
 
    hash_stock_name = {}
    hash_name.each do |llave, valor|
      hash_id.each do |llave1, valor1|
        if llave == valor1 
          hash_stock_name[llave1] = valor
        end
      end
    end
  
    total = (@f2_wild.to_f) + (@f2_p1p2.to_f) + (@f2_p1.to_f) + (@f2_p2.to_f)
    wild_expected = (9.0/16) * total
    p1_expected = (3.0/16) * total
    p2_expected = (3.0/16) * total
    p1p2_expected = (1.0/16) * total
    wild_observed = (@f2_wild.to_f)
    p1_observed = (@f2_p1.to_f)
    p2_observed = (@f2_p2.to_f)
    p1p2_observed = (@f2_p1p2.to_f)
    @chi_2 = (((wild_observed - wild_expected)** 2) / wild_expected) + (((p1_observed - p1_expected)** 2) / p1_expected) + (((p2_observed - p2_expected)** 2) / p2_expected )+ (((p1p2_observed - p1p2_expected)** 2) / p1p2_expected)
    
    #If the Chi-square value calculated for an experiment is greater than that corresponding to the 5% probability, the hypothesis is rejected.
    #The hypothesis is that genes are not linked
    #5% for 3 degrees is 7.82
    if @chi_2.to_f >= 7.82
      @linked = "These genes are linked"
      hash_stock_name.each do |llave2, valor2|
        if @parent1 == llave2
          @gene_name1 = valor2
        elsif @parent2 == llave2
          @gene_name2 = valor2
        end
      end
      #Report messages 
      puts "Recording: #{@gene_name1} is genetically linked to #{@gene_name2} with chisquare #{@chi_2}"
      puts "Final Report:\n#{@gene_name1} is linked to #{@gene_name2}\n#{@gene_name2} is linked to #{@gene_name1}"
      linked_array << @gene_name1
      linked_array << @gene_name2
    end
  end 
end