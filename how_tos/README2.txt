#Running simc and passing an input file automatically

#Store the input file name in a file "input"
	echo current.data > input
	./simc < input

#Alternatively do:
	./simc <<EOF 
	current.data
	EOF
