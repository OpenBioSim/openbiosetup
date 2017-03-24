#!/usr/bin/python
#
#
# Driver script for biomolecular simulations structures setup pipeline
#
#

class Sequence():
    """
    Describes the sequence of a molecule to include in the system
    being prepared.

    Sequence is an abstract base class. Specialised classes implement
    'organic', 'protein', 'dna', 'rna' sequences.
    """

    def __init__:
        pass


class Structure():
    """Describes the coordinates, topology, residue and atom identifiers present
    in a molecule structure.

    Structure is an abstract base class. Specialised classes are
    used for 'organic', 'protein', 'dna', 'rna' structures.

    """
    def __init__:
        pass

    def coordinates:
        """Returns the coordinates of the structure"""
        pass

    def connectivity:
        """Returns the bond/bond-orders of the structure"""

    def sequence:
        """Returns the sequence described by the structure"""
        pass

class Model():
    """
    Describe the sequence and structure of a molecule to include in
    the system being prepared.

    Model is an abstract base class. Specialised classes are
    used for 'organic', 'protein', 'dna', 'rna' models.

    Remark:
      clarify differences between Model and Structure classes

    """
    def __init__:
        pass

def loadSequences([]):
    """
    Input:
    A list of strings that contains paths to input files

    Action:
    Load sequences contained in list of input files,

    Code has to handle various sequence input formats
    Initially we can rely on file extensions.
    We store different types of sequences in different categories
    e.g. first fasta file read -->
                      if it contains a protein, store in proteinsequence1
                      if it contains dna, store in dnasequence1
          first smiles file read --> organicsequence1
    We should also do some sanity checks (is this a valid SMILES etc..?)

    classify sequence type and create sequence objects

    Allowed sequence types are

    organicX --> described by SMILES
    proteinX --> described by IUPAC/IAB one-letter amino acid code
    dnaX     --> described by IUPAC/IAB one-letter base-pair code
    rnaX     --> described by IUPAC/IAB one-letter base-pair code

    Output:

    A sequences dictionnary e.g.

    sequences = { "protein1" : ProteinSequenceObject ,
    "organic1": OrganicSequenceObject }"""

    return {}

def loadStructures([]):
    """
    Input:
    A list of strings that contains paths to input files

    We must support a variety of file formats (PDB, mol2, ...)

    Action:

    Load structures contained in list of input files.
    Classify structure types and create structure objects

    Remark:
      * structural waters, ions .... handled as 'organic' structures

    Output:

    A structures dictionnary e.g.

    structures ={ "protein1" : ProteinStructureObject,
                  "protein2" : ProteinStructureObject,
                  "organic1" : OrganicStructureObject }

    """

    return {}

def mapProteinSequences(sequences, structures):
    """Input:
    A system dictionnary

    Action:
    Find all passed protein sequences
    For each protein sequence
    - create a new proteinmodel variable (ProteinStructureClass)
    perform pairwise sequence alignment to sequences present in all
    protein structures passed

    If exact match (100% identity)
          use coordinates of protein structure
    if sequence is fully contained within template structure (substructure):
          use coordinates of protein structure
    if partial match (similarity > minimum threshold)
         pick template protein sequence with highest similarity
         use homology modelling to construct structures of query sequence
    if no match or similarity < minimum threshold
         BLAST sequence to find possible templates (need access to database)

    Remark:
       if no match, it may be challenging to model quarternary structures so
       should consider bailing out if more than one protein sequence needs a
       structure for this setup

       Some protein structures may have no matching sequences. Default behavior
       is to ignore structures.

    Output:
    A list of proteinmodelX entries

    each proteinmodelX variable is a ProteinModelobject
    """
    return []

def reviewProteinModels([]):
    """
    Input:
    A list of protein models

    Action:
    For each model, create a scene to render the structure of the model
    Allow user to edit structure coordinates


    Output:
    A list of protein models

    """
    return []



def mapOrganicSequences(sequences, structures):
    """
    Input:
    A list of sequences and structures objects

    Action:
    For each organic sequence in the passed sequences
    create an organicmodel object

    Do a fingerprints similarity calculation against all
    sequences in organic structures

    full match (Tanimoto 1):
      use coordinates of matched sequence for model coordinates.

    partial match (Tanimoto > threshold):
       perform unsupervised alignment onto coordinates of matched sequence.

    no match:
       generate 3D coordinates from SMILES, activate flag for random
       position/orientation upon embedding

    Output:
    A list of organicmodelX entries.
    each organicmodelX variable is a OrganicModelObject
    """
    return []

def reviewOrganicModels([]):
    """
    Input:
    A list of organic models

    Action:
    For each model, create a scene to render the structure of the model
    Allow user to edit structure coordinates


    Output:
    A list of organic models

    """
    return []

def embedModels(proteinmodels, organicmodels):
    """
    Input:
    A list of proteinmodel and organicmodel objects

    Action:
    Create an empty 'space' container
    Add every protein and organic models that have known
    'absolute' coordinates to that container

    Define system bounding box

    Compute bounding box of each positioned protein and organic
    objects

    For each remaining model lacking absolute coordinates
    Compute bounding box of model

    Keep randomly inserting model in a random orientation
    into system bounding box, until no clashes between molecules
    bounding boxes.

    Remark:
       Allow flexibility in bounding box geometry

    Output:
    A system dictionnary that contains a list of protein and organic models,
    and information about the bounding box geometry
    """
    return {}

def solvate(system, environment):
    """
    Input:
    A system dictionnary that contains models within a bounding box
    An environment dictionary that defines the (co)solvents

    Action:
    Work out number of each cosolvent molecule to insert in system bounding box
    to achieve concentrations specified by environment

    Here need alg to position cosolvent molecules without significant clashes

    Remark:
      - We get more realistic distributions by replicating
        pre-equilibrated solvent boxes

    Output:
    A system dictionnary with organic models for each cosolvent molecule
    """
    return system

def adjustProtons(system, pH=7.0):
    """
    Input:
    A system dictionary containing models of molecules that
    have titrable sites
    A pH value

    Action:
    Algorithm that determines protonation states from structures

    Output:
    A system dictionary with updated titrable sites

    """
    return system


def writeOutput(system):
    """Input:
    A system dictionnary

    Action:
    Write the contents of the system to output file formats.

    We need to output
    - Molecules (connectivity, coordinates, bond orders, elements/atom types)
    - Box dimensions

    """


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
    ###################################################################
    #
    # The following steps are core to the pipeline
    #
    # We load the input sequences
    sequences = loadSequences(sequenceFiles)

    # We load the input structures
    structures = loadStructures(structureFiles)

    # Depending on the input the next steps may not be always executed.
    # Multiple different implementations of each action may become available.
    # Suggests need a different design than top-level module function for
    # the pipeline 'Actions'

    # We map protein sequences onto protein structures to make protein models
    proteinmodels = mapProteinSequences(sequences, structures)
    # * separate procedures could be implemented for DNA/RNA sequences

    # Optional human interaction to review generated models, select
    # one conformation among ambiguous solutions, or manually edit
    # structures
    # commented out as not needed for minimal implementation
    # proteinmodels = ReviewProteinModels(proteinmodels)

    # We map organic sequences onto organic structures
    organicmodels = mapOrganicSequences(sequences, structures)

    # Optional human interaction to review generated models, select
    # one conformation among ambiguous solutions, or manually edit
    # structures
    # commented out as not needed for minimal implementation
    # organicmodels = ReviewOrganicModels(organicmodels)

    # We embedd all proteinmodels and organicmodels in a system
    system = embedModels(proteinmodels, organicmodels)

    # We solvate the system
    system = solvate(system, environment)

    # We optimise protonation states
    # commented out as not needed for minimal implementation
    # system = adjustProtons(system, pH=environment['pH'])

    # Optional human interaction to review ambiguous protonation states

    #
    # The following step is core to the pipeline
    #
    # We write this information in files suitable as input
    # to a simulation engine input prep utility
    writeOutput(system)
    ###################################################################
    ###################### Pipeline ends ##############################
    ###################################################################
