Setup:

Upgrade gocmds:

`sudo gocmd upgrade`
`gocmd init`

You will be prompted to insert your CyVerse credentials and then you can run:

`gocmd get --progress /iplant/home/shared/NCEMS/working-groups/oceans-of-disorder/minimal-database/minimal_noenv.db.gz .`

And finally unpack the database: `gunzip minimal_noenv.db.gz`

Now, open up the jupyter notebook and follow the instructions to make a query based on a COG ID, plot and compare the results by IDR property, and save the extracted data to a .csv file (if desired)
