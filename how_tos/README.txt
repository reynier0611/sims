Get latest simc version from Werner
	git clone https://github.com/boeglinw/simc_gfortran.git

Go into the downloaded directory and "make"
	cd simc_gfortran/
	make clean;make

Go into the root directory and "make"
	cd root/
	make
	make -f fmake_tree.mak

Create a directory outside the simc main directory
and setup simc there
	cd ../..
	cd ..
	mkdir temp_dir

Copy the setup script from the simc main directory and
past it into the new work directory
	cp simc_gfortran/setup.sh /temp_dir
	cd temp_dir

Edit the setup script, and give it the path to the
work directory
	vim setup.sh 

Source the setup script
	sh setup.sh 

That should create all necessary links. Now we can
run simc from the work directoty
	./simc

