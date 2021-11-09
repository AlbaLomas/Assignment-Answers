require 'rest-client'
require 'json'
require './Interaction_Network.rb'

#Enter the input file (gene list)
gene_text = ARGV[0]

#Enter the threshold input (a number)
threshold = ARGV[1]

#Have the list of genes in gene list
@lista_gene = []
file=File.open(gene_text,"r")
file.each do |line|
	line = line.strip
	gene = line.downcase
	@lista_gene << gene
end


#Count number of genes of gene list
file2=File.open(gene_text,"r")
a = file2.readlines.size

#######################################################################################################################################################
#Initilize Interaction_Network.rb file

Interaction_Network.open_genelist(gene_text)
Interaction_Network.search_interactors(threshold)
Interaction_Network.createdirect_network()
Interaction_Network.createindirect_network()
Interaction_Network.modificate_indirect_network()
Interaction_Network.join_nets()
Interaction_Network.get_kegg_go()
#Interaction_Network.show_results()
Interaction_Network.create_class()
Interaction_Network.numberofnetworks
Interaction_Network.get_net

##############################################################################################################################################################################

#Know which genes in the different networks are listed genes
included_genes = Array.new
Interaction_Network.class_variable_get(:@@interactionnetwork_array).each do |a|
	a.members_list.each do |g|
		#Not add to included_genes list if it is already in 
		if @lista_gene.include? (g) 
			included_genes << g 

		end
	end
end

#Remove repeated genes again
included_genes = included_genes.uniq
#puts included_genes
#Num of genes involved in indirect and direct networks
num_ingen = included_genes.length()
#to float
num_ingen = num_ingen.to_f
#total genes of arabidopsis file
total_genes = a.to_f


#Percentage of genes in the list that interact with each other
@total = (num_ingen/total_genes) * 100

open("report.txt", "w") do |archivo|

	archivo.puts("REPORT FILE\n")
	archivo.puts("\n")
	archivo.puts("This document reports which list members of ArabidopsisSubNetwork_GeneList.txt interact with one another, in networks of two or three genes with the following structure: A - C / A - B - C, A and C are listed genes and B is a common interactor.\n")
	archivo.puts("\n")
	archivo.puts("\t - It would only be considered that 'A' interacts with 'B' and vice versa when the intact-miscore is equal or greater than the threshold that is indicated by the user. It is recommended to use a threshold of 0.485 (the one proposed in the following paper --> DOI: 10.1371/journal.pone.0108567).\n")
	archivo.puts("\n")
	archivo.puts("This document also reports the KEGG and GO functional annotations of those interacting members.\n")
	archivo.puts("\t- It should be emphasized that it would be annotated both KEGG ID and Pathway Name.\n")
	archivo.puts("\t- It should be emphasized that it would be only annotated the GO Terms associated with the biological process part (GO:ID and GO Term Name).\n ")
	archivo.puts("\n")
	#Number of analysed genes
	archivo.puts("The number of genes that are contained in ArabidopsisSubNetwork_GeneList.txt, and therefore are being analysed is #{a}.\n")
	#Number of networks
	archivo.puts("The analysis shows that there are #{Interaction_Network.numberofnetworks} direct and indirect networks of genes.\n")
	archivo.puts("\n")
	#Write the members of the different networks and the interactor gene
	@contador_direc_net = []
	@contador_ind_net = []
	Interaction_Network.class_variable_get(:@@interactionnetwork_array).each do |red| #vamos iterando todas las networks obtenidas
		archivo.puts("Network number #{red.network_num}.\n")
		archivo.puts("\tThe genes that conform this interaction network are: ")
		archivo.puts(red.members_list)
		###########################################################################################################################################
 	
 		#Indirect networks 
	 	if red.members_list.length() == 3
	 	
	 		archivo.puts("This is an indirect network")
	 		archivo.puts("\n")
	 		red.members_list.each do |member|
	 			#Count the number of genes involved in direct networks 	
				if @lista_gene.include? (member)
					@contador_ind_net << member
				else
					@interactor = member
					archivo.puts("\tIn the network number #{red.network_num} the interactor gene is #{@interactor}.\n")
					archivo.puts("\n")
				end
			end
			@contador_ind_net = @contador_ind_net.uniq
			@contador_ind = @contador_ind_net.length()

		#Direct networks 

	 	else
	 		red.members_list.each do |member2|
	 			#Count the number of genes involved in direct networks 	
				if @lista_gene.include? (member2)
					@contador_direc_net << member2
				end
			end
			@contador_direc_net = @contador_direc_net.uniq
			@contador_direc = @contador_direc_net.length()
			
			archivo.puts("\n")
	 		archivo.puts("This is a direct network")
	 	end
	 	#########################################################################################################################################
	 	
	 	archivo.puts("\n")
	 	archivo.puts("KEGG ANNOTATION.\n")
		archivo.puts("\n")

		#Write the annotation if it exists
		if red.kegg_pathways.empty? 
			archivo.puts("None of the genes that integrate the network #{red.network_num} have a KEGG Pathways annotation associated with them.\n")
			archivo.puts("\n")
		else
			archivo.puts("The research shows that the genes of network #{red.network_num} have #{red.kegg_pathways.length()} KEGG Pathways annotations.\n")
			archivo.puts("\n")
				
		#Not all network genes have an annotation
			red.kegg_pathways.each do |annotation_k|
				archivo.puts("\t-KEGG ID: #{annotation_k[0]} / - KEGG Pathway: #{annotation_k[1]}\n")
				archivo.puts("\n")
					
			end

		end

		#GO
			
		archivo.puts("GO ANNOTATION.\n")
		archivo.puts("\n")

		#Write the annotation if it exists
		if red.go_terms.empty? 
			archivo.puts("None of the genes that integrate the network #{red.network_num} have GO Terms associated.\n")
			archivo.puts("\n")
		else
			archivo.puts("The research shows that the genes of network #{red.network_num} have #{red.go_terms.length()} GO Terms associated.\n")
			archivo.puts("\n")

		#Not all network genes have an annotation
				
			red.go_terms.each do |annotation_go|
				archivo.puts("\t-GO ID: #{annotation_go[0]} / -GO Biological process: #{annotation_go[1]}\n")
				archivo.puts("\n")

			end
		end
	end
	archivo.puts("\n")
	total_ind = (@contador_ind.to_f / total_genes) * 100
	total_dir = (@contador_direc.to_f / total_genes) * 100
	archivo.puts "The percentage of listed genes that are involved in indirect networks is #{total_ind}"
	archivo.puts("\n")
	archivo.puts "The percentage of listed genes that are involved in direct networks is #{total_dir}"
	archivo.puts("\n")
	archivo.puts "The percentage of listed genes that are involved in different networks is #{@total}"

	
end
