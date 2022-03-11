#!/usr/bin/env python3

"""
Title: pathchecker.py
Created on Thur 09-03-2022
Author: Daniel Nilsson

Description:
    Utilities for handling path.
    * Creation of directories
    * check if a file already exists


List of functions:
    - check_path_overwrite   
    - createDirs
    
List of "non standard" modules:


"""

from pathlib import Path

#Path essentials


def check_path_overwrite(path):
    """
    Function checking if given output path exists and if it is a file. 
    If already present ask the user if it should be overwritten
    Returns True if file should be overwritten. Otherwise False
    """
    excecute_script = True #variable that holds a boolean for overwriting decision
    path = Path(path) #turn the input argument to a Path object
    
    #check if path does exist and if it is a file
    if path.exists() and path.is_file():
         
         #If file exists. Ask for premission to overwrite the file
          usr_input_overwrite = input("\nSpecified output file does already exist." + "\n" +
                                 "Do you want to overwrite it? (y/n): ")
        
          #Based on user answer, set excecute_script to True or False
          if usr_input_overwrite == "y": #if y is chosen the set 
              excecute_script = True     #excecute_script to True
                   
          elif usr_input_overwrite == "n": #if n is chosen the set 
              excecute_script = False      #excecute_script to False
                
          else: #if something else than y or n then set excecute_script to False
              excecute_script = False 
              
    #If the file does not exist, return True       
    else:
       excecute_script = True
       
    return excecute_script #Return a boolean (True or False)




def createDirs(filename):
    """
    This first function only does one thing. That is to strip the file name of the path that should be created.
    If this were to be done in the function loopDirs(), there would be a problem with removing a path each recursion.
    This would cause directories to be missed. This is solved by wrapping the inner function by an additional function so the command 
    for removing the file name is only done once.   
    """

    def loopDirs(path):
        """ 
        This function will take a Path instance as an input and then checks
        if it exists. If not, save the given non-existing directory name 
        and traverse up one level to check if the parent directory
        exists and so on until an existing directory path is found.
        Then check if it is an absolute path. If not, join current working directory
        with the relative path.
        Lastly it creates the new path by creating each directory at a time  
        """
        
        dirs_to_create = []  #a list to hold the directories to create
        
        #check if directory exists   
        if not path.exists(): #if not existing 
           dirs_to_create.append(path.name) #append last directory to the list  
           path = loopDirs(path.parent) #call itself to check if next parent exists 
           
           if not path.is_absolute(): #Checking if the given path is absolute
               path = Path.cwd().joinpath(path) #append path to current working dir if not absolute
                                                #to get a complete path to write
            
                     
        for path_to_make in dirs_to_create: #iterate over the list of directory names to create
           path = path.joinpath(path_to_make) #join path to create with existing path
           path.mkdir() #create the new directory
           
           
        return path #returns new output path
    
    loopDirs(Path(filename).parent) #strip the file name before calling the loopDirs() to create the path (directories) to the output file
    


    
if __name__ == "__main__":
    """
    
    If the script is run as main.
    Here is the number of arguments provided checked and if they are valid.
    If everything looks fine then the functions are called  
    
    """ 
    pass