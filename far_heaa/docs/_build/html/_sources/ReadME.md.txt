# Fast and Robust High Entropy Alloy Analysis

## Features

- Use DFT generated mixing enthalpy values to create regular solution models that can be used to generate phase diagrams.
- A visualization toolkit that plots binary, ternary phase diagrams along with other higher order visualizations.
- High throughput calculations across phase space for stable High Entropy Alloys to assist experimentalists with design choices.

## Installation

```commandline
git clone
```

```commandline
cd src/far_heaa
pip install -r requirements.txt
```
## Setup

 - ### Step 1
    Open the database/metadata.json and peek inside for tags that need to be given to the project. <br>
   (add metadata description)
 - ### Step 2
    The default metadata needs only your materials project api-key to be used. If you don't have the api-key as yet, switch the im_flag to false and run the examples provided down. 
 - ### Step 3
    Peek Inside the implementations folder. All the use cases have been provided. It is recommended to create a new folder in the farheaa and write your own code with loops, conditions and other custom code to chain these implementations as required. <br>
    ##### Note: Once you run the implementation, the visualization will be stored in a newly created plots folder, with a unique directory structure.
 - ### Step 4
    If your implementations require changing the metadata flags frequently, then there are handlers to update the metadata and change the flags.
````commandline
from far_heaa.io.metadata_handler import MetadataHandler

mH = MetadataHandler()

# print the metadata keys
mH.access_metadata_keys()

# update a certain key
mH.update_metadata(key: 'end_member', value: 'xyz')

meta_data = mH.get_metadata
````

For more details, there is always documentation. 
## Change Log

## Cite

## Usage Policy

