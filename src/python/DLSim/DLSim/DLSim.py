import os
from types import SimpleNamespace
import pandas as pd

# load user defined functions
try:
    from func_lib import (load_cpp_shared_library,
                          check_required_files_exist,
                          get_filenames_from_folder_by_type,
                          func_running_time,
                          path2linux
                          )

except Exception:
    from .func_lib import (load_cpp_shared_library,
                           check_required_files_exist,
                           get_filenames_from_folder_by_type,
                           func_running_time,
                           path2linux
                           )


__all__ = ['DLSim']


class DLSim:
    """Perform network assignment using the selected assignment mode

     Parameters
        assignment_mode
            0: Link-based UE, only generates link performance file without agent path file

            1: Path-based UE, generates link performance file and agent path file

            2: UE + dynamic traffic assignment (DTA), generates link performance file and agent path file

            3: ODME

        column_gen_num
            number of iterations to be performed before on generating column pool

        column_update_iter
            number of iterations to be performed on optimizing column pool

    """

    def __init__(self,
                 assignment_mode: int = 2,
                 column_gen_num: int = 10,
                 column_update_num: int = 10) -> None:
        self.assignment_mode = assignment_mode
        self.column_gen_num = column_gen_num
        self.column_update_num = column_update_num

        # define current directory as working directory
        self.DLSim_WORKING_DIR = path2linux(os.getcwd())

        # define the required files
        self.DLSim_required_inputs = ["settings.csv", "node.csv", "link.csv", "demand.csv"]

    def check_working_directory(self, working_dir: str = os.getcwd()) -> None:
        """Check the working directory.

        Args:
            path_dir (str, optional): The path of the working directory. Defaults to os.getcwd().

        Returns:
            None: Change the python working environment to a specified working_dir.

        Examples:
            >>> from DLSim import DLSim as DL
            >>> DL.check_working_directory()
                Current working directory: C:\\Users\\luoxiangyong\\Desktop\\DLSim

            >>> DL.check_working_directory("C:\\Users\\luoxiangyong\\Desktop\\DLSim\\test")
                Current working directory: C:\\Users\\luoxiangyong\\Desktop\\DLSim\\test
        """
        # prepare stardardard cross-platform path
        working_dir = path2linux(working_dir)

        # change the working directory to the user specified working directory
        os.chdir(working_dir)

        # assign the user specified working directory to the global variable
        self.DLSim_WORKING_DIR = working_dir

        # remind the user to put all the required files in the working directory
        print(f"Current working directory: {self.DLSim_WORKING_DIR}. \
            Please note that all the input files {self.DLSim_required_inputs} should be in this directory. \
            You can change the working directory by calling the function: \
            check_working_directory(working_dir = 'your working directory')")

    def check_DLSim_input_files(self):
        """Check the required input files in the working directory.

        Examples:
        >>> from DLSim import DLSim as DL
        >>> DL.check_DLSim_input_files()
            ALL the required files ['settings.csv', 'nodes.csv', 'link.csv', 'demand.csv'] in the working directory.

        # if the required files are not in the working directory
        >>> DL.check_DLSim_input_files()
            Please prepare missing files and put them in the working directory
        """

        files_in_working_dir = get_filenames_from_folder_by_type(self.DLSim_WORKING_DIR, "csv")

        isRequiredFilesExist = check_required_files_exist(self.DLSim_required_inputs, files_in_working_dir)

        if isRequiredFilesExist:
            print(f"All the required files {self.DLSim_required_inputs} in the working directory.")
        else:
            print(f"Please prepare missing files and put them in the working directory: {self.DLSim_WORKING_DIR}.")

    @property
    def DLSim_settings(self) -> SimpleNamespace:
        """Load the settings from the settings.csv file.

        # The current logic is to open the settings.csv file first,
        # load all parameters,
        # and update parameter if necessary,
        # then, write the updated parameters back to the settings.csv file

        # The reason to do so is the cpp shared library read the settings.csv file only.
        # we can not pass parameters to the cpp shared library directly at current stage.


        Returns:
            SimpleNamespace: a SimpleNamespace object with all the parameters in the settings.csv file. (TODO)

        Examples:
            >>> from DLSim import DLSim as DL
            >>> DL.check_working_directory()
            >>> DL.check_DLSim_input_files()
            >>> DL.DLSim_settings
        """

        path_settings = path2linux(os.path.join(self.DLSim_WORKING_DIR, "settings.csv"))

        print(f"Loading settings from {path_settings}. \
            you can prepare or change the settings.csv file in the working directory.")

        return pd.read_csv(path_settings)

    @func_running_time
    def perform_kernel_network_assignment_simulation(self):
        """Perform kernel network assignment simulation.

        Examples:
            >>> from DLSim import DLSim as DL
            >>> DL.check_working_directory()
            >>> DL.check_DLSim_input_files()
            >>> DL.perform_kernel_network_assignment_simulation()
        """

        print("The default assignment mode is 2 (UE + DTA) \
               and the default column generation number is 10 \
               and the default column update number is 10.")

        print("assignment_mode: ", self.assignment_mode)
        print("column_gen_num: ", self.column_gen_num)
        print("column_update_num: ", self.column_update_num)

        # make sure assignment_mode is right
        assert self.assignment_mode in {0, 1, 2, 3}
        # make sure iteration numbers are both non-negative
        assert self.column_gen_num >= 0
        assert self.column_update_num >= 0

        print('\nDTALite run starts')

        # load the cpp shared library
        _dtalite_engine = load_cpp_shared_library()

        _dtalite_engine.network_assignment(self.assignment_mode,
                                           self.column_gen_num,
                                           self.column_update_num)

        print('\nDTALite run completes\n')
        print(
            f'check link_performance.csv in {os.getcwd()} for link performance\n'
            f'check agent.csv in {os.getcwd()} for unique agent paths\n'
        )

if __name__ == "__main__":

    # load the DLSim class
    DL = DLSim()

    path_test = r"C:\Users\roche\Anaconda_workspace\001_Github\DLSim-MRM\datasets\ASU"

    # check the working directory
    DL.check_working_directory()

    # check all the required files exist
    DL.check_DLSim_input_files()

    # load and update settings
    DL.DLSim_settings

    # user can change setting parameters here
    DL.setting.time_period = "0700_0800"

    # perform kernel network assignment simulation
    DL.perform_kernel_network_assignment_simulation()
