#Create an “InteractionNetwork” Object to contain the members of each network
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
	@@genes = []
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

	#Create an array for storaging all the features of the Interaction_Network class
	@@interactionnetwork_array = []

	#List with all the genes in the MARK file that we will extract with the def open_genelist
	
  
	def initialize(params={}) 
		@network_num = params.fetch(:network_num, 0) 
		@members_list = params.fetch(:member_list, [])
		@kegg_pathways = params.fetch(:kegg_pathways, [])
		@go_terms  = params.fetch(:go_terms, [])
		@@interactionnetwork_array << self #Introduce every attribute of every class 
		@@number_networks += 1 #Count the number of networks
	end

  	   
  def self.fetch(url, headers = {accept: "*/*"}, user = "", pass="")
  #class method to retrieve web data avoiding possible crashes
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
 
  
	  
	#Function to list all the genes in the in @@genes
	def self.open_genelist(genelistfile) #el self es necesario o no funciona!!!
		File.open(genelistfile).each_line do |line|
			line = line.strip #to no create lists
      line = line.downcase
			@@genes << line
		end
	
	end
	
	def self.search_interactors(threshold=0.485, genes=@@genes) #class method to search the interactions of the genes from the previous file
    #and save them in a class variable if MIscore is bigger than the cutoff value (threshold)
    #New hash for the interacting relations, key (listed gene) and the values are lists of the listed gene interactors 
    	@@interacting=Hash.new

    	genes.each do |gene|

        #Change the gene name to catch the genes in tab 25
    		gene2 = "tair:" + gene.downcase
    		#Hash for the interactors of listed gene
        lista = []

    	 #search each listed gene in utoronte database
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
  
            next if score2.to_f <= threshold #Only catch those interactions with a score equal or greater than the threshold 
          			
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
      #contador = 0
      #@@interacting.each do |key, values|
        #contador +=1 
        #puts key
        #puts values
      #end
    #puts contador 
  end  
  def self.createdirect_network()
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
            netdir = netdir.sort
            if !@@list_dirgen.include? (netdir)
              @@list_dirgen << netdir
            else
              next
            end
          end
        end
      end
    end

    @@list_dirgen = @@list_dirgen.uniq
    number1 = @@list_dirgen.length()
    puts "we have #{number1} networks of two genes"
    
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
    @@list_indgen = @@list_indgen.uniq
    number = @@list_indgen.length()

		puts "we have #{number} networks of three genes"


  end


  def self.join_nets()
    @@nets = @@list_indgen + @@list_dirgen
    number4 = @@nets.length()
    puts number4
  end
   

    ###Hasta aquí parece que va bien


  def self.get_kegg_go()
    #parameters(db is the database )
    #Method for annotate Kegg Pathways 
    #Array of arrays with all the Kegg Pathways wich are associated with the different networks
    
    @@dic_final = Hash.new
    @@nets.each do |net|
      lista_kegg = []
      lista_go = []
      net.each do |gene|
        web = Interaction_Network.fetch("http://togows.dbcls.jp/entry/kegg-genes/ath:#{gene}/pathways.json")
        if web
          ans = JSON.parse(web.body)[0] #access the results (json format is a list)
          ans.each do |kegg_id,kegg_p| #data is in hash format, get
                lista_kegg << [kegg_id,kegg_p] #save the annotations of each locus in the list
          end
        end
        #get go 
        resp = Interaction_Network.fetch("http://togows.org/entry/uniprot/#{gene}/dr.json")
        if resp
          ans2=JSON.parse(resp.body)[0] #definimos el diccionario
          ans2["GO"].each do |go| 
            spl = go[1].split(":")
            if spl[0] == 'P'
              go_id = go[0].sub 'GO:', ""
              go_p = go[1].sub 'P:', ""
              go_2 = [go_id, go_p]
              if lista_go.include? go_2
                next
              else
                lista_go << go_2  
              end
            end
          end
        end
      end 
      @@dic_final[net] = [lista_kegg,lista_go]
        #return lista_kegg #return the list of annotations for the network
   
    end
  end

  def self.show_results()
    
    @@dic_final.each do |key, values|
      puts key
      puts values[0]
      puts values[1]
    end
  end

  def self.create_class
    n_net=0 #for the network number
    @@dic_final.each do |net2, annotations| #interacting class variable is a hash
      n_net+=1 #add one each time a new class is created 
      Interaction_Network.new(:network_num => n_net, :members_list =>net2, :kegg_pathways =>annotations[0], :go_terms =>annotations[1])
    end
  end

  def self.numberofnetworks #get the total number of networks
    return @@number_networks  
  end

  def self.get_net # get the list with all the objects
    return @@interactionnetwork_array
  end

end



  