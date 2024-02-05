# Bayesian-Analysis-of-Shopping-Habits-with-Gamma-Poison-Model
This project aims to provide insights into the shopping habits of customers in a small business by utilizing Bayesian analysis and the Gamma-Poisson model.

#Steps to run the project.

1. Open RStudio.
2. Create a new project:
    a. Select "File" > "New Project".
    b. Choose "Existing Directory".
    c. Navigate to the directory where your R script and data file are located.
    d. Click "Create Project".
3. In the RStudio environment, navigate to the "Files" pane in the lower right corner of the window.
4. Open the R script file:
    a. Double-click on it.
    b. Select it and click "Open".
5. To run the code, install the below libraries beforehand.
    a. install.packages("MASS")
    b. install.packages("MCMCpack")
  Note: Gelman-Rubin diagnostic has been run both with and without the package. Calculations have been done without explicitly using
        the predefined package.
6. Data have been directly fed into the code itself as our final data is derived from the original dataset.
7. Took the original dataset and derived secondary data from it which will act as our priors. so, no need for original data to import.
8. Click "Run" on the right corner to run the code.
