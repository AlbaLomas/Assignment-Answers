require "date"
#class seed-stock

class Seedstock 
  # create "attribute accessor"
  @@all = []
  attr_accessor :seed_stock 
  attr_accessor :mutant_gene_id
  attr_accessor :last_planted
  attr_accessor :storage
  attr_accessor :grams_remaining 
  
  
  def initialize (params = {}) # get a name from the "new" call, or set a default
    @seed_stock = params.fetch(:seed_stock, '0000')
    @mutant_gene_id= params.fetch(:mutant_gene_id, "AT0G00000")
    @last_planted = params.fetch(:last_planted, "00/00/0000")
    @storage = params.fetch(:storage, "cama0")
    @grams_remaining = params.fetch(:grams_remaining, "0")
  end

  #Method for actulize the remaining amount of seed after planting an amount (number) and date which would be the day where the file is updated
  def actualize_rem_seed_and_date (number)
    @grams_remaining = @grams_remaining - number
    @last_planted = DateTime.now.strftime('%d/%m/%Y')
    #The amount of seed can not be less than zero, if it is, it will appear a warning message
    if @grams_remaining <= 0
       @grams_remaining = 0
       puts "WARNING: we have run out of Seed Stock #{@seed_stock}"
    end
    
  end

  def test_array()
    puts "The seed has #{@seed_stock} and there are #{grams_remaining} remaining grams"
  end
  
  #Method for writing the actualized file 
  def write_file()
    "#{@seed_stock}\t#{mutant_gene_id}\t#{last_planted}\t#{@storage}\t#{grams_remaining}\n"
  end
end