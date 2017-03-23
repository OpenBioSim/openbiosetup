#!/usr/bin/python
#
#
# Driver script for biomolecular simulations structures setup pipeline
#
#

class Structure():
    """Placeholder to store information about
    coordinates, topology, residue and atom identifiers present
    in a molecule
    """
    def __init__:
        pass


def loadSequences([]):
    """Load sequences contained in list of input files,
    and classify sequence type.

    Allowed sequence types are

    organicX --> described by SMILES
    proteinX --> described by IUPAC/IAB one-letter amino acid code
    dnaX     --> described by IUPAC/IAB one-letter base-pair code
    rnaX     --> described by IUPAC/IAB one-letter base-pair code

    Output:

    A sequences dictionnary e.g.

    sequences = { "protein1" : "NAC..." ,
    "organic1": "CC(O)C" }"""

    return {}

def loadStructures([]):
    """Load structures contained in list of input files,
    and classify structure types

    Output:

    A structures dictionnary e.g.

    structures ={ "protein1" : structure_object,
                  "protein2" : structure_object,
                  "organic1" : structure_object }

    """

    return {}

if __name__ == '__main__':

    ###################################################################
    # User input data
    # We assume the user knows which sequences & molecules
    # they want to structural models for !
    # Here we pass fasta files to describe protein sequences and
    # canonical smiles for organic molecules
    # The implicit assumption is that one sequence = one molecule
    sequenceFiles = ["humanthrombin-lightchain.fasta",
                     "humanthrombin-heavychain.fasta",
                     "3J.smi"]
    # The user also specifies the environment that will surround
    # the structures derived from the input sequences
    # Arbitrarily complex buffers are described by adding
    # multipled cosolvent entries
    # We use disconnected smiles to group together components
    # into electroneutral groups
    # * For colsolvents that contain titrable residues (e.g. Tris below)
    # we will have to come up with suitable heuristics to generate
    # electroneutral boxes
    # One drawback of describing water as a cosolvent is that we will not
    # have geometries matching exactly a given rigid body water model that
    # we may wish to use later on for a simulation.
    # We are going to say for now that this is a problem for rigid body water
    # forcefields, and deal with that later down the pipeline
    # Note setting a pH value is meaningful for aqueous solutions only
    environment = {'pH':7.4,
                   'cosolvent1':('O','55 M')
                   'cosolvent2':('[Na+].[Cl-]','50 mM'),
                   'cosolvent3':('OCC(N)(CO)CO','20 mM') }

    # The user also specifies any structural data that should be used
    # by the pipeline to generate models for the input sequences.
    # Here we load a single PDB file
    structureFiles = ["2ZC9.pdb"]
    ###################################################################
    ###################### Pipeline begins ############################
    # We initially define our system as a flat python dictionary
    # Later on when we understand what type of data we need to add & save
    # we will explore other file formats
    # Initially the system is empty
    system = {}
    # We add to the system the environment information
    for key, value in environment.items():
        system[key] = value

    # We load the sequences
    # Code has to handle various sequence input formats
    # Initially we can rely on file extensions.
    # We store different types of sequences in different categories
    # e.g. first fasta file read -->
    #                  if it contains a protein, store in proteinsequence1
    #                  if it contains dna, store in dnasequence1
    #      first smiles file read --> organicsequence1
    # We should also do some sanity checks (is this a valid SMILES etc..?)
    sequences =  loadSequences(sequenceFiles)
    # All sequences are stored in the system
    for sequenceid, content in sequences.items():
        system[sequenceid] = content

    # We load the input structures. We must support a variety of input
    # format.
    # Code returns a dictionnary of structures, each entry corresponds to
    # one molecule.
    # For our purposes at this stage a structure is a collection of
    # atomic coordinates and bond orders. This will be kept in a class.
    structures = loadStructures(structureFiles)
    # All structures are stored in the system
    for structureid, content in structures.items():
        system[structureid] = content

    # Next steps
    # output of mapping is a 'model' e.g. proteinmodel1, organicmodel1 ..
    # * map protein sequences onto protein structures
    # work out how to deal with partial and no matches
    # partial matches are resolved by homology modelling using partial match as
    # template
    # warn if partial match is so bad that might as well not use input
    # structure?
    # no matches can also be handled via fully automated homology modelling
    #  (if there are no matches though we can't predict quaternary structures
    #   current code may limit itself to single sequences with no match)
    #
    # * separate procedure to be implemented for DNA/RNA sequences (not now)
    # * map organic sequences onto organic structures
    #   full match: use structure coordinates
    #   partial match: perform unsupervised alignment onto
    #   most chemically similar match
    #   no match: guess 3D coordinates from SMILES and randomly insert in
    #   system (but not uniformly, bias to avoid clashes) This requires box
    #   volume definition. Might be useful to generate N random insertions
    #   for some setups (e.g. MSM binding )
    #
    # * solvate systems
    #   work out bounding box that gives satisfactory cutoff, insert cosolvents
    #
    # * optimise protons
    #   best done before or after cosolvents placement ?
    # 
