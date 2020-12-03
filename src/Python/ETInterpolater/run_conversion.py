from ETInterpolater import ETInterpolater, ETQuantities
from mpi4py import MPI
from time import sleep
import argparse
import datetime
import os
from shutil import copy


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
            "geometry_resolution": 400,
            "nz": 1,
            "nr": [0],
            "nt": [0],
            "np": [0],
            "r_limits": [0.0],
            "it": [0],
            "pickle": 0,
            "pickle_folder": "",
            "interpolation": 2,
            "quantities": ["alp"]
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
                if len(words) != 2:
                    raise ValueError("The line %s is not correctly formatted" %line)
                else:
                    read_parameters[words[0].strip()] = words[1].strip()

        self._evaluate_parameters(read_parameters)



    def _evaluate_parameters(self, read_parameters):
        for key in read_parameters.keys():
            if key in self.parameters.keys():

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
        if len(self.parameters["r_limits"]) != self.parameters["nz"]:
            raise ValueError("r_limits does not have a length (%s) of nz+1 (%s)"
                %(len(self.parameters["r_limits"]), self.parameters["nz"]+1))


    



class Setup:
    def __init__(self, file_name):
        self.p_object = ParameterReader(file_name)
        self.parameters = self.p_object.parameters
        self.datetime = datetime.datetime.now()
        self.parameters["time"] = str(self.datetime)
        self.folder = self.parameters["result_folder"] + "/results_%s" %self.datetime.strftime("%d_%m_%y_%H.%M.%S")
        
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
            copy(from_c_path +".C", self.folder+"/")
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
            f.write(str(self.parameters["nz"]) + "\n")
            for n in ["r_limits","nr", "nt", "np"]:
                for i in self.parameters[n]:
                    f.write(str(i) + " ")
                
                if n != "np":
                    f.write("\n")





class Runner:
    def __init__(self, args):
        self.args = args
        if args.option == "run":
            self._run()

    def _run(self):
        parameters_loc = self.args.parameters
        if os.path.isfile(parameters_loc):
            setup = Setup(parameters_loc)
            setup.create_setup()
        else:
            raise ValueError("Parameter file not found at: %s" %parameters_loc)

        folder = setup.parameters["simulation_folder"]
        pickle_folder = setup.parameters["pickle_folder"]
        pickle = setup.parameters["pickle"]
        quantites = setup.parameters["quantities"]
        it = setup.parameters["it"]
        geometry_corner = setup.parameters["geometry_size"]
        geometry_res = setup.parameters["geometry_resolution"]
        #c_path = setup.parameters["c_path"]
        result_path = setup.parameters["result_path"]
        c_path = result_folder

        inter = ETInterpolater(folder, setup.parameters["interpolation"])
        g = inter.make_positive_geometry([-geometry_corner]*3, geometry_res)
        et_q = ETQuantities(g, it, folder, pickle_folder=pickle_folder, pickle=pickle)
        
        inter.analyse_bbh(g, et_q, it, result_path=result_path,c_path=c_path, quantities=quantites, test=False)
    
    def _continue(self):
        pass
    
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
    parser.add_argument("-f", "--folder", dest="restult_folder", help="Folder where old results are found, in case continue, finish or rerun is choosen")
    parser.add_argument("-p", "--parameters", dest="parameters", default="parameters.txt", help="Name of parameter file for run. If nothing is given, it is assumed to be in the same location as the python program, and named parameters.txt")


    args = parser.parse_args()
    runner = Runner(args)


    #conversion_setup = Setup("parameters.txt")
    #conversion_setup.create_setup()
    