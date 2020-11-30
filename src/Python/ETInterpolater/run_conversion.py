from ETInterpolater import ETInterpolater, ETQuantities, ParameterReader
from mpi4py import MPI
from time import sleep
import argparse
import datetime
import os
from shutil import copy


"""
TODO:
- Create argpasers so the user can start new or continue old conversion, or go final gyoto_conversion
- Create function to check if more quantities have been added to already existing conversion
- Make safty checks if dz, nr, nt or np have been changed when continuing conversion
- Make argparses to do all the things 
    - Full/partial conversion
    - Final gyoto conversion, with check if all quantites are there
    - test plots with python
    - test plots with lorene?
- Make a function that gets the origin from ET, both so user can see it
    -> Function for final gyoto conversion.
"""



class Setup:
    def __init__(self, file_name):
        self.p_object = ParameterReader(file_name)
        self.parameters = self.p_object.parameters
        self.datetime = datetime.datetime.now()
        self.parameters["time"] = str(self.datetime)
        self.folder ="./" + self.parameters["result_folder"] + "/results_%s" %self.datetime.strftime("%d_%m_%y_%H.%M.%S")
        print(self.parameters["result_folder"])
        
    def create_setup(self):
        if not os.path.isdir(self.folder):
            os.mkdir(self.folder)
        else:
            raise ValueError("Folder already exists")

        self.create_parameterfile()

        from_c_path = ""

        if self.parameters["c_path"] == "":
            from_c_path = "./get_points"
        elif "get_points" in self.parameters["c_path"][-1]:
            from_c_path = self.parameters["c_path"]
        else:
            from_c_path = "./" + self.parameters["c_path"] + "/get_points"

        try:
            copy(from_c_path, self.folder+"/")
        except:
            raise ValueError("get_points not found at: %s" %from_c_path)
        
    
    def create_parameterfile(self):
        with open(self.folder + "/parameterfile.txt", "w") as f:
            for key in self.parameters.keys():
                if self.parameters[key] == "":
                    f.write(key + ": " + '""' + "\n")
                else:
                    f.write(key + ": " + str(self.parameters[key]) + "\n")
        
        with open(self.folder + "/lorene_parameters.txt", "w") as f:
            f.write(str(self.parameters["dz"]) + "\n")
            for n in ["nr", "nt", "np"]:
                for i in self.parameters[n]:
                    f.write(str(i) + " ")
                
                if n != "np":
                    f.write("\n")


    




    
        








class Runner:
    def __init__(self):
        pass




if __name__ == "__main__":
    conversion_setup = Setup("parameters.txt")
    conversion_setup.create_setup()