scripts
=======

master's thesis 

This is a MATLAB class library for sturcturing hippocampal place field data. It also includes several methods to perform 
unit analysis and others to aid visualization. The purpose was to develop a class library to support different datasets.
Currently, the "kenji" and "MTA" datasets are supported. To add more, the @GenericTrial/Convert2Generic.m has to be updated. 

-------

Conventions:

        variables - start with lower case letters
                    exception, logical variables have all upper case letters
        methods   - start with upper case letters
    
    
    
The data has to be located in the ~/data/<dataset>/ .

All results are saved to the ~/<dataset>/analysis folder.

For MTA dataset, the original tree structure is adopted for compatibility.

-------

Class names:
  
        GenericTrial - The base class. This instances of this class have properties and methods to aid unit analysis. 
        GenericPF    - Place field calss. This class has additional properties and methods to assess place fields. 
        
-------

Usage:

        trial = GenericTrial('filebase', 'trialname');
        trial = GenericTrial('filebase');
        
        the trialname is prompted if not specified, all the trialnames are displayed for convenience
      
        Example : gt = GenericTrial('ec013.931_942', 'ec013.932')
        
                     
------
        
For additional help, check " doc GenericTrial " or " doc GenericPF ". 

