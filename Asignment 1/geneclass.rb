#class gene_information

class Geneinformation
  @@gene_all = []
  attr_accessor :gene_id
  attr_accessor :gene_name
  attr_accessor :mutant_phenotype
  attr_accessor :linked
  
  
  def initialize (params = {}) # get a name from the "new" call, or set a default
    @gene_id = params.fetch(:gene_id, "AT0G00000")
   	#lo de some gene se puede cambiar  
    @gene_name= params.fetch(:gene_name, "some gene")
    @mutant_phenotype = params.fetch(:mutant_phenotype, [])
    @linked = params.fetch(:linked, "No linked")
  end

  def new_atributes()
    @@gene_all << self
    @@gene_all << gene_id
    @@gene_all << gene_name
    @@gene_all << mutant_phenotype
  end
  def linked_(ufo,pi)
    if ufo .to_s == @gene_name.to_s
      @linked == "#{@gene_name} is linked to #{pi}"
    elsif pi .to_s == @gene_name.to_s
      @linked == "#{@gene_name} is linked to #{ufo}"
    end
  
  end 
  def test_()
  puts "#{@gene_id}\t#{@gene_name}\t#{@mutant_phenotype}\t#{@linked}"
  end
end