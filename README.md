# EIC-Dipper  

Documentation:  
https://www.overleaf.com/read/fxvswvktwmfd
https://www.overleaf.com/read/rtjrybdvgprj#db7d9b

# How to install  

1. Clone this repository:  
```
$ git clone https://github.com/yhoffmann/EIC-Dipper.git
```
2. Enter the directory and clone the submodules:
```
$ cd EIC-Dipper/
EIC-Dipper$ git submodule update --init --recursive
```

# How to build
1. Create a directory for the obj binaries
```
EIC-Dipper$ mkdir obj
```
2. Build the submodules, libs and externals (this will probably only need to be done once per pull):
```
EIC-Dipper$ make libs
```
3. Build the main part of the code (also use this every time you want to build after modifying something that is not a submodule):
```
EIC-Dipper$ make
```
