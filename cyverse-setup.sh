echo "Upgrading gocmd..."
sudo gocmd upgrade
echo "Initializing gocmd..."
gocmd init
echo "Getting the database..."
gocmd get --progress /iplant/home/shared/NCEMS/working-groups/oceans-of-disorder/minimal-database/minimal_noenv.db.gz .
echo "Unpacking the database..."
gunzip minimal_noenv.db.gz
echo "DONE."
