# Query the ocean IDR database

* This GitHub repository enables easy searching of the SQLite database `minimal_noenv.db` using a Jupyter notebook.

## Getting started

1. Kicking off a CyVerse session

* Navigate to `de.cyverse.org` on your preferred web browser
* On the navigation menu on the left, click the `Apps` tab
* Click on `JupyterLab Datascience` to start the process of launching a session
![Apps menu page](images/Apps-menu.png)
* You will now need to navigate a series of screens to choose your parameters and launch the session:
* (1) "Analysis Info" - on this screen, just click the `Next` button in the bottom righthand corner
![Analysis Info page](images/analysis-info.png)
* (2) "Advanced Settings" - here we need to select the computational resources for the job. Select `CPU Cores = 4`, `Minimum Memory = 16.0 GiB`, and `Minimum Disk Space = 64 GB` using the dropdown menus. 
![Advanced Settings page](images/advanced-settings.png)
* (3) "Review and Launch" - click the `Launch Analysis` button in the bottom righthand corner of the screen
![Review and Launch](images/review-and-launch.png)
* (4) On the next screen, hit `Go to Analysis` at the top
![Go to Analysis](images/go-to-analysis.png)
* You should now see a loading bar like the one below. It may take a few minutes to get your resources set up. If it stalls, try refreshing the page. 
![Launching the app](images/launching-app.png)

2. Clone the repository

Note well, you will need to set up a GitHub access token. If you have not done so already, see the tutorial **here**. 

3. Configuring `gocmd` and getting the database file

Setup:

Upgrade gocmds:

`sudo gocmd upgrade`
`gocmd init`

You will be prompted to insert your CyVerse credentials and then you can run:

`gocmd get --progress /iplant/home/shared/NCEMS/working-groups/oceans-of-disorder/minimal-database/minimal_noenv.db.gz .`

And finally unpack the database: `gunzip minimal_noenv.db.gz`


3. Making a query


Now, open up the jupyter notebook and follow the instructions to make a query based on a COG ID, plot and compare the results by IDR property, and save the extracted data to a .csv file (if desired)

