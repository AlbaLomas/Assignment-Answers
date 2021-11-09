#Create an “Interaction_Network” Object to contain the members of each network
#Annotate it with any KEGG Pathways the interaction network members are part of
#	both KEGG ID and Pathway Name
# Annotate it with the GO Terms associated with the total of all genes in the network
#	BUT ONLY FROM THE biological_process part of the GO Ontology!
# 	Both GO:ID and GO Term Name

# The objective is checking the interactions between the listed Arabidopsis genes.

class Interaction_Network 
	#First, record/count the number of networks (starts in 0)
	#With each object this number increases 1
	@@number_networks = 0 

  #different list that we will need after
	@@genes = [] #list of al genes of the Arabidopsis list 
  @@list_gen = []
  @@listat = []
  @@nets = []

	#Then, create a variable which contains the number of each network
	attr_accessor :network_num
	
	#Create a variable that contains a list with the genes inside each network
	attr_accessor :members_list

	#Annotate any KEGG Pathways, this variable is a list for each network with the KEGG pathways of the members list
	attr_accessor :kegg_pathways

	#Annotate the GO Terms associated with the total of all genes in the network
	attr_accessor :go_terms

	#Create an array for storaging all the features of the Interaction_Network class for the report
	@@interactionnetwork_array = []


	
  
	def initialize(params={}) 
		@network_num = params.fetch(:network_num, 0) #An integer
		@members_list = params.fetch(:members_list, 'NA') # A list
		@kegg_pathways = params.fetch(:kegg_pathways, 'NA') # A list 
		@go_terms  = params.fetch(:go_terms, 'NA') # A list
		@@interactionnetwork_array << self #Introduce every attribute of every class when a class is created 
		@@number_networks += 1 #Count the number of networks
	end

  	   
  def self.fetch(url, headers = {accept: "*/*"}, user = "", pass="")
  #Retrive webs and avoid error or problems 
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
      return response  
    rescue RestClient::Exception => e
      $stderr.puts e.inspect
      response = false
      return response  
    rescue Exception => e
      $stderr.puts e.inspect
      response = false
      return response  
  end
 
  
	  
	#Function to asign all the genes to the @@genes list
	def self.open_genelist(genelistfile) # self is necessary or the script not works!!
		File.open(genelistfile).each_line do |line|
			line = line.strip #to no create lists
      line = line.downcase #(downcase) to match with the results of the web page
			@@genes << line 
		end
	
	end
	
	def self.search_interactors(threshold, genes=@@genes, specie = "taxid:3702") 
    #class method to search the interactions of the genes from the previous file
    #and save them in a class variable if MIscore is bigger than the cutoff value (threshold) which is fixed by the user
    #New hash for the interacting relations, key (listed gene) and the values are lists of the listed gene interactors 
    	@@interacting=Hash.new

    	genes.each do |gene|

        #Change the gene name to catch the genes in tab 25
    		gene2 = "tair:" + gene.downcase
    		#Hash for the interactors of listed gene
        lista = []

    	 #search each listed gene in utoronto database
    		page = Interaction_Network.fetch("http://bar.utoronto.ca:9090/psicquic/webservices/current/search/query/#{gene}/?format=tab25")

        #if there is a response
   			if page
   				#Obtain the different tables (Tab25 Format) by \n splitting 
   				page = page.split("\n")
   				page.each do |a|

            a = a.split("\t")
            score = a[14] #threshold in the 15th position 
            gene_name = a[2].downcase #gene in the 3th or 4th position
            genename2 = a[3].downcase
            
            score2 = score.sub "intact-miscore:", "" #Only catch the threshold number to filter
  
            next if score2.to_f <= threshold.to_f #Only catch those interactions with a score equal or greater than the threshold 
          	next if !a[9] == specie || !a[10] == specie #Only catch those interactions if the genes belong to the specie
          	#Find the listed gene and its interactor and catch only the gene name 
            if gene2 == gene_name

              g_h = a[2].sub "tair:", ""
						  g_2 = a[3].sub "tair:", ""
              g_2 = g_2.downcase
					
						  lista << g_2

            elsif gene2 == gene_name && gene2 == genename2
              next

            else
						  g_h = a[3].sub "tair:", ""
						  g_2 = a[2].sub "tair:", ""
              g_2 = g_2.downcase

						  lista << g_2
					
					
   					end
   		
   				end
          		#If there list is not empty store it in the hash
   				if lista.empty?  == false
   				
   					@@interacting[gene.downcase] = lista	
    			end
   			
    		end

    	end

      # Test the @@interacting hash
      #contador = 0
      #@@interacting.each do |key, values|
        #contador +=1 
        #puts key
        #puts values
      #end
    #puts contador 
  end  

  def self.createdirect_network() #Catch direct networks between listed genes

    @@list_dirgen = []

    #Go through the interacting hash
    #This loop will search the interactors of each listed gene and if there is a common interactor between two listed genes the network will be created
    @@interacting.each do |key, values|
      values.each do |value|
        if key == value
          next
        else
          if @@genes.include?(value)
            netdir = [key,value]
            #Sort the genes of the network due to delete repeated networks 
            netdir = netdir.sort
            #Not include repeated networks
            if !@@list_dirgen.include? (netdir)
              @@list_dirgen << netdir
            else
              next
            end
          end
        end
      end
    end

    #Make sure that there are not repeated networks
    @@list_dirgen = @@list_dirgen.uniq
    number1 = @@list_dirgen.length()

    #puts "we have #{number1} networks of two genes"
    
    #@@list_dirgen.each do |net|
      #puts net
    #end
    
  end

  def self.createindirect_network()
    #This method will create the A-B-C networks between the gene in the interacting hash, A and C are listed genes. 

    #Create a list which would be replaced with the different networks of three genes
    @@list_indgen = []

    #Go through the interacting hash
    #This loop will search the interactors of each listed gene and if there is a common interactor between two listed genes the network will be created
    @@interacting.each do |key, values|
    	llave1 = key
    	values.each do |value|
    		value2 = value
    		@@interacting.each do |key2, values2|
    			if key2 == llave1
    				next
          elsif key2 == value
            next

    		  else
    				values2.each do |value3|
    					value4 = value3
              if value2 == value4
                cosass = [llave1, value2, key2]
                #Sort the genes of the network due to delete repeated networks 
                cosass = cosass.sort
                #If the list is repeated it would not be included
                if @@list_indgen.include? (cosass)
                  next

                else
                  @@list_indgen << cosass
                end
    				  end
            end
    		  end
        end
      end
    end

    #Make sure that there are not repeated networks
    @@list_indgen = @@list_indgen.uniq

    #Know the number of indirect networks
    #number = @@list_indgen.length()
    #puts "we have #{number} networks of three genes"


  end

  def self.modificate_indirect_network() 
    #It is possible that two listed genes have different interactors so there would be redundant networks. 
    #We have to search if the different list have 2 listed genes in common and delete this redundant networks. 
    #If the repeated genes are not listed, the network would not be deleted as one of the repeted genes is an interactor so the networks are different. 
  
    @@list_indgen.each do |list|
      #list.each do |g|
      @@list_indgen.each do |list2|
        next if list == list2
        if list2.include? (list[0]) 
          if list2.include?(list[1]) 
            if @@genes.include? (list[0]) and @@genes.include? (list[1])
              @@list_indgen.delete(list)
            else
              next
            end

          elsif list2.include?(list[2]) 
            if @@genes.include? (list[0]) and @@genes.include? (list[2])
              @@list_indgen.delete(list)
            else
              next
            end

          else 
            next

          end

        elsif list2.include? (list[1]) 
          if list2.include?(list[0])
            if @@genes.include? (list[1]) and @@genes.include? (list[0])
              @@list_indgen.delete(list)
            else
              next
            end

          elsif list2.include?(list[2]) 
            if @@genes.include? (list[1]) and @@genes.include? (list[2])
              @@list_indgen.delete(list) 
            else
             next
            end

          else
            next
          end

        elsif list2.include? (list[2]) 
          if list2.include?(list[0]) 
            if @@genes.include? (list[2]) and @@genes.include? (list[0])
              @@list_indgen.delete(list)
            else
              next
            end
          
          elsif list2.include? (list[1]) 
            if @@genes.include? (list[2]) and @@genes.include? (list[1])
              @@list_indgen.delete(list)
            else
              next
            end
          end
        end
      end
    end
    #Test the definitive number of networks of three genes
    #number7 = @@list_indgen.length()   
    #puts "we have #{number7} networks of three genes"
    #@@list_indgen.each do |genes|
      #puts genes
    #end
  end
         
  def self.join_nets() #Join the direct and indirect networks to search kegg and go terms
    @@nets = @@list_indgen + @@list_dirgen
    number4 = @@nets.length()
    #Know the total number of networks
    #puts number4
  end
   
  def self.get_kegg_go()

    #Method for annotate Kegg Pathways and go terms
    #@@n_net will return the number of each network and will be the keys of the @@dict_final hash, the values of each key will be a list of lists in which the first position corresponds with the network members; the second, with the kegg pathways; and the third, with the go terms of each gene of the network. 
    @@n_net=0
    @@dic_final = Hash.new
    @@nets.each do |net|
      lista_g = [] #Create an empty list for each network to store the genes that belong to the network and add it to the dict_final hash
      @@n_net+=1  #With each networks n_net increases 1
      lista_kegg = [] #Create an empty list for each network to store the kegg pathways of all genes of the network
      lista_go = [] #Create an empty list for each network to store the go annotations of all genes of the network

      net.each do |gene|
        lista_g << gene
        #Search each gene of the networks in kegg web 
        web = Interaction_Network.fetch("http://togows.dbcls.jp/entry/kegg-genes/ath:#{gene}/pathways.json")
        if web
          ans = JSON.parse(web.body)[0] #json format hash
          ans.each do |k_id,k_p| #data is in hash format
                lista_kegg << [k_id,k_p] 
          end
        end

        #get go 
        #Search each gene of the networks in go web 
        resp = Interaction_Network.fetch("http://togows.org/entry/uniprot/#{gene}/dr.json")
        if resp
          ans2=JSON.parse(resp.body)[0] 
          ans2["GO"].each do |go| #json format hash
            #Only catch GO Terms associated with the biological process part (GO:ID and GO Term Name)
            spl = go[1].split(":")
            if spl[0] == 'P'
              go_id = go[0].sub 'GO:', ""
              go_p = go[1].sub 'P:', ""
              go_2 = [go_id, go_p]

              #Do not include if the Go term is repeated
              if lista_go.include? go_2
                next
              else
                lista_go << go_2  
              end
            end
          end
        end

      # Fill the final dict
      @@dic_final[@@n_net] = [lista_g,lista_kegg,lista_go]
      end 
        #return lista_kegg #return the list of annotations for the network
   
    end
  end

  def self.show_results() #Show the hash that would be insert in the classes, only used to test
    
    @@dic_final.each do |key, values|
      puts key
      #puts values[0]
      #puts values[1]
    end
  end

  def self.create_class #method to create the classes
    @@dic_final.each do |key2, annotations| #interacting class variable is a hash
     Interaction_Network.new(:network_num => key2, :members_list =>annotations[0], :kegg_pathways =>annotations[1], :go_terms =>annotations[2])
      
    end

  end

  def self.numberofnetworks #get the total number of networks
    return @@number_networks  
  end

  def self.get_net(all = @@interactionnetwork_array) # get the list with all the objects
    return all
  end
end



  