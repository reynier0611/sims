******************************************
              ---  SIMC  ---
       Notes by Reynier Cruz Torres
******************************************

I got this simc version from Werner Boeglin (boeglinw@fiu.edu)
It works slightly different from the default simc


******************************
Steps to run this simc version
******************************
In the infiles directory, define your input file and name it "current.inp"

From the main directory do:
> make clean
> make
> ./simc

This should ask for the input file mentioned above. Thus, type:

> current.inp

simc will start running


******************************
Handling the output file
******************************
When simc finishes, the output will be an ascii file named "current.data"
From the main directory run the python script make_root_one_default.py 
in the python directory:

> python python/make_root_one_default.py

This code has the name "current.data" hard-coded in. That's why above
I said the input file had to be renamed to current.inp

This code converts the ascii file into a root file

******************************
Making an ntuple reading code
******************************
In the same directory where "current.root" is located do the following:

> root -l
> f = new TFile("current.root")
> SNT->MakeClass()

SNT stands for "Simc N-Tuple"
A .C and a .h files are created. These can be used to get information from
the ntuple. This last "Making an ntuple reading code" is a default thing in
root. I just added it here as a reminder for myself
