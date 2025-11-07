## Table of Contents

1. [Installation](#installation)
2. [Project Structure](#project-structure)
3. [Running the Code](#running-the-code)
4. [License](#license)

## Installation

This code runs based on **MATLAB 2022b**. Ensure that you have this version installed on your system.

1. **Clone the repository**:
   - To ensure stable downloading even with an unstable network connection, you can directly go to the URL **[https://github.com/leiyuk/O2UDI](https://github.com/leiyuk/O2UDI)**, click **'Download ZIP'**, and then unzip the file.
2. **Navigate to the project folder**:
    - Open MATLAB on your computer.
    - In MATLAB, navigate to the folder where the O2UDI project is located.


## Project Structure

The project is divided into three main folders:

### Folder 1: `simulation`
- **Description**: The simulation validation for the evolution of debris clouds from explosive breakup, corresponding to the first part of the Results section in the paper.
- **Files**:
    - `main.m`: The main program, which performs the debris cloud evolution calculations using parallel computing and stores the results in `.mat` files.
    - `plot_NC.m`: This script generates plots that display the comparison between the reference values and inferred values, as well as the JS divergence.
    - `plot_gif_new.m`: This script outputs a video that dynamically shows the comparison between the reference values and inferred values throughout the process.
    - Other files: Custom functions called during the execution of the program.


### Folder 2: `single_debris`
- **Description**: Validates the core support of the O2UDI method using TLE data from the Cosmos 1408 breakup event. Corresponds to the second part of the Results section in the paper.
- **Files**:
    - `main.m`: The main program, which outputs the linear fitting results of the semi-major axis change rate and B* based on data from the Cosmos 1408 breakup event.
    - `cosmos_1408`: Contains the TLE data of 100 debris pieces generated from the Cosmos 1408 breakup. All debris data was obtained from the website [https://www.space-track.org/](https://www.space-track.org/).
    - Other files: Custom functions called during the execution of the program.






### Folder 3: `multi_debris`
- **Description**: Validates the O2UDI method using TLE data from the Cosmos 2251 breakup event. Corresponds to the third part of the Results section in the paper.
- **Files**:
    - `main.m`: The main program that, after reading 10 years of debris cloud data, infers the spatial distribution of small debris based on large-scale debris in the cloud, and compares it with the true distribution.
    - `plot_NC.m`: This script generates plots comparing the true values with the inferred values of the distribution, as well as the calculated JS divergence.
    - `plot_gif.m`: This script outputs a video that dynamically shows the comparison between the true values and inferred values throughout the process.
    - `cosmos_2251`: Contains TLE data of 288 debris pieces generated from the Cosmos 2251 breakup, from the time of breakup to 10 years later.
    - Other files: Custom functions called during the execution of the program.



## Running the Code

### Folder 1: `simulation`
1. Open MATLAB and navigate to the `simulation` directory. Then, run the `main.m` file. Due to the large number of debris, the computation is slow. For example, using 76 cores of an Intel 8368 CPU, the calculation takes about 3 hours.  

2. Run `plot_NC.m` and `plot_gif_new.m` to generate the figures and video.


### Folder 2: `single_debris`
1. Open MATLAB and navigate to the `single_debris` directory. Run the `main.m` file to generate the linear fitting results, and the relative errors of R² and the intercept will be displayed in the command window.



### Folder 3: `multi_debris`
1. Open MATLAB and navigate to the `multi_debris` directory. Run the `main.m` file. The calculation takes about 10 minutes on an 8-core CPU computer, with most of the time spent reading the TLE data.  

2. Run `plot_NC.m` and `plot_gif.m` to generate the figures and video.

## License

This repository is released under the **Apache License 2.0**.  
Users are permitted to use, reproduce, and modify the code for research and educational purposes, provided that proper attribution is given to the original authors.



To see the full license text, visit [http://www.apache.org/licenses/LICENSE-2.0](http://www.apache.org/licenses/LICENSE-2.0).
