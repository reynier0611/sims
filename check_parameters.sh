rm test_file.data
ln -s $1 test_file.data
./test_parameters < $1
