require 'rest-client'
require 'json'
require './prueba3.rb'

Interaction_Network.open_genelist("text.txt")
Interaction_Network.search_interactors()
Interaction_Network.createdirect_network()
Interaction_Network.createindirect_network()
Interaction_Network.get_kegg_go()
#Interaction_Network.show_results()
Interaction_Network.create_class
Interaction_Network.numberofnetworks
Interaction_Network.load
Interaction_Network.get_net.each do |a|
	b = a.members_list
	puts b
end
	

