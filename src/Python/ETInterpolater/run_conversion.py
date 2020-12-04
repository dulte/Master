from ETInterpolater import ETInterpolater, ETQuantities
from mpi4py import MPI
from time import sleep
import argparse
import datetime
import os
from shutil import copy
import numpy as np


"""
TODO:
- [~] Create argpasers so the user can start new or continue old conversion, or go final gyoto_conversion
- Make so analyze_bbh used the paths from the parameter file
- Create function to check if more quantities have been added to already existing conversion
- Make safty checks if nz, nr, nt or np have been changed when continuing conversion
- Make argparses to do all the things 
    - Full/partial conversion
    - Final gyoto conversion, with check if all quantites are there
    - test plots with python
    - test plots with lorene?
- Make a function that gets the origin from ET, both so user can see it
    -> Function for final gyoto conversion.


- Implement non positive geometry

- Make the code use all the iterations...
"""

class ParameterReader:
    
    def __init__(self, file_name):
        self.parameters = {
            "time": 0,
            "c_path": "",
            "simulation_folder": "",
            "result_folder": "",
            "geometry_type": "positiv",
            "geometry_size": 20,
            "geometry_resolution": 50,
            "nz": 1,
            "nr": [0],
            "nt": [0],
            "np": [0],
            "r_limits": [0.0],
            "it": [0],
            "pickle": 0,
            "pickle_folder": "",
            "interpolation": 2,
            "quantities": ["alp"],
            "finish": 0,
        }
        self.file_name = file_name
        self._read_file()


    def _read_file(self):
        read_parameters = {}
        with open(self.file_name) as f:
            for line in f:
                if line == "\n":
                    continue
                words = line.split(":")
                if len(words) != 2 and words[0].strip() != "time":
                    raise ValueError("The line %s is not correctly formatted" %line)
                else:
                    read_parameters[words[0].strip()] = words[1].strip()

        self._evaluate_parameters(read_parameters)



    def _evaluate_parameters(self, read_parameters):
        for key in read_parameters.keys():
            if key in self.parameters.keys():
                if key == "time":
                    continue

                if type(eval(read_parameters[key])) == type(self.parameters[key]):
                    self.parameters[key] = eval(read_parameters[key])
                else:
                    raise ValueError("For the parameter %s, the type %s did not match the expected %s" 
                        %(key, type(eval(read_parameters[key])), type(self.parameters[key])))

        self._give_warnings(read_parameters)    
        

    def _give_warnings(self, read_parameters):

        for p in self.parameters.keys():
            if p not in read_parameters.keys():
                if p in ["np", "nt", "nr", "nz", "r_limits"]:
                    raise ValueError("%s is not found. This must be given" %p)
                elif "folder" in p:
                    self.parameters[p] = "./"                
                elif p == "time":
                    continue
                else:
                    print "[~] %s not given, running with default value %s"%(p,self.parameters[p])

        if self.parameters["pickle_folder"] == "" and self.parameters["pickle"]:
            _ = raw_input("[~] Pickle is True, but no folder is given for pickling. This means that the pickles are saved here. Press any key to continue...")
        

        
        if len(self.parameters["nr"]) != self.parameters["nz"]:
            raise ValueError("nr does not have a length (%s) of nz (%s)"
                %(len(self.parameters["nr"]), self.parameters["nz"]))
        if len(self.parameters["np"]) != self.parameters["nz"]:
            raise ValueError("np does not have a length (%s) of nz (%s)"
                %(len(self.parameters["np"]), self.parameters["nz"]))
        if len(self.parameters["nt"]) != self.parameters["nz"]:
            raise ValueError("nt does not have a length (%s) of nz (%s)"
                %(len(self.parameters["nt"]), self.parameters["nz"]))
        if len(self.parameters["r_limits"]) != self.parameters["nz"]+1:
            raise ValueError("r_limits does not have a length (%s) of nz+1 (%s)"
                %(len(self.parameters["r_limits"]), self.parameters["nz"]+1))


    



class Setup:
    def __init__(self, file_name):
        print "[~] Reading Parameters"
        self.p_object = ParameterReader(file_name)
        self.parameters = self.p_object.parameters
        self.datetime = datetime.datetime.now()
        self.parameters["time"] = str(self.datetime)
        self.folder = self.parameters["result_folder"] + "/results_%s" %self.datetime.strftime("%d_%m_%y_%H.%M.%S")
        print "[+] Reading Parameters Done"
        
    def create_setup(self):
        print "[~] Setting up enviroment"
        if not os.path.isdir(self.folder):
            os.mkdir(self.folder)
        else:
            raise ValueError("Folder already exists")

        self.create_parameterfile()
        self._make_positionfile()

        from_c_path = ""

        if self.parameters["c_path"] == "":
            from_c_path = "./get_points"
        elif "get_points" in self.parameters["c_path"][-1]:
            from_c_path = self.parameters["c_path"]
        else:
            from_c_path = "./" + self.parameters["c_path"] + "/get_points"

        try:
            copy(from_c_path, self.folder+"/")
            #copy(from_c_path +".C", self.folder+"/")
        except:
            raise ValueError("get_points not found at: %s" %from_c_path)

        os.chdir(self.folder)
        self.folder = "./"

        print "[+] Setting up enviroment done"
        
    
    def create_parameterfile(self):
        with open(self.folder + "/parameterfile.txt", "w") as f:
            for key in self.parameters.keys():
                if self.parameters[key] == "":
                    f.write(key + ": " + '""' + "\n")
                elif type(self.parameters[key]) == str:
                    f.write(key + ": " + '"' + str(self.parameters[key]) + '"' + "\n")
                else:
                    f.write(key + ": " + str(self.parameters[key]) + "\n")
        
        with open(self.folder + "/lorene_parameters.txt", "w") as f:
            f.write(str(self.parameters["nz"]) + "\n")
            for n in ["nr", "nt", "np", "r_limits"]:
                for i in self.parameters[n]:
                    f.write(str(i) + " ")
                
                if n != "r_limits":
                    f.write("\n")
    
    def load_setup(self, folder):
        #self.parameters["c_path"] = "./"
        #self.parameters["result_folder"] = "./"
        self.folder = folder
        os.chdir(self.folder)
        self.folder = "./"
        self._read_positionfile()




    def _make_positionfile(self):
        inter = ETInterpolater(self.parameters["simulation_folder"], self.parameters["interpolation"])
        possible_iterations, positions, radii = inter.read_bbh_diag()
        with open(self.folder + "/bh_data.txt", "w") as f:
            f.write(str(radii[0,0]) + "\n")
            f.write(str(radii[1,0]) + "\n")
            for i in positions[0,:,0]:
                f.write(str(i)+", ")
            f.write("\n")
            for i in positions[1,:,0]:
                f.write(str(i)+", ")

    def _read_positionfile(self):
        with open(self.folder+"/bh_data.txt", "r") as f:
            r = f.readline().split()[0]
            r2 = f.readline().split()[0]
            self.radii = [float(r), float(r2)]
            self.pos1 = []
            self.pos2 = []
            for i in f.readline().split(",")[:-1]:
                self.pos1.append(float(i))
            for i in f.readline().split(",")[:-1]:
                self.pos2.append(float(i))
            



class Runner:
    def __init__(self, args):
        self.args = args
        if args.option == "run":
            self._run()
        elif args.option == "continue":
            self._continue()

    def _find_setup(self):
        pass

    def _setup(self):
        pass





    def _run(self):
        parameters_loc = self.args.parameters
        if os.path.isfile(parameters_loc):
            setup = Setup(parameters_loc)
            setup.create_setup()
        else:
            raise ValueError("Parameter file not found at: %s" %parameters_loc)

        exit()
        folder = setup.parameters["simulation_folder"]
        pickle_folder = setup.parameters["pickle_folder"]
        pickle = setup.parameters["pickle"]
        quantites = setup.parameters["quantities"]
        it = setup.parameters["it"]
        geometry_corner = setup.parameters["geometry_size"]
        geometry_res = setup.parameters["geometry_resolution"]
        c_path = setup.folder + "/"
        result_path = setup.folder
        do_gyoto_converstion = setup.parameters["finish"]
        

        inter = ETInterpolater(folder, setup.parameters["interpolation"])
        g = inter.make_positive_geometry([-geometry_corner]*3, geometry_res)
        et_q = ETQuantities(g, it[0], folder, pickle_folder=pickle_folder, pickle=pickle)
        
        inter.analyse_bbh(g, et_q, it, result_path=result_path,c_path=c_path, quantities=quantites, test=False, do_gyoto_converstion=do_gyoto_converstion)
    
    def _continue(self):
        result_folder = self.args.result_folder
        if not os.path.isdir(result_folder):
            raise ValueError("%s directory not found" %result_folder)
        
        setup = Setup(result_folder+"/parameterfile.txt")

        setup.load_setup(result_folder)

    
    def _rerun(self):
        pass

    def _testplot(self):
        pass
    
    def _finish(self):
        pass




if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument("option", choices=["run", "continue", "rerun", "testplot", "finish"], help="Choose action: run, continue, rerun, finish or testplot")
    parser.add_argument("-q", "--quantity", dest="plot_quantity", help="Quantity to be plotted if testplot is choosen")
    parser.add_argument("-f", "--folder", dest="result_folder", help="Folder where old results are found, in case continue, finish or rerun is choosen")
    parser.add_argument("-p", "--parameters", dest="parameters", default="parameters.txt", help="Name of parameter file for run. If nothing is given, it is assumed to be in the same location as the python program, and named parameters.txt")


    args = parser.parse_args()
    runner = Runner(args)


    #conversion_setup = Setup("parameters.txt")
    #conversion_setup.create_setup()
    