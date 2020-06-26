# Make errors fatal
set -e

# Package to retrieve Google Drive files
pip install gdown

# Download and unzip the data
cd $TRAVIS_BUILD_DIR
gdown https://drive.google.com/uc?id=$DRIVE_FILE_ID
unzip -o data.zip
