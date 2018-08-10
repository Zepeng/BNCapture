$Id: README 80190 2014-04-07 10:18:04Z gcosmo $

     =========================================================
     Geant4 - an Object-Oriented Toolkit for Simulation in HEP
     =========================================================

                            Hadr03
                            ------

   How to compute total cross section from the direct evaluation of the 
   mean free path ( see below, item Physics).
   How to identify nuclear reactions.
   How to plot energy spectrum of secondary particles.	 
	
 1- GEOMETRY DEFINITION
 
   It is a single box representing a 'semi infinite' homogeneous medium.
   Two parameters define the geometry :
 	- the material of the box,
	- the (full) size of the box.
 	
   The default geometry (10 m of molybdenum) is built in DetectorConstruction,
   but the above parameters can be changed interactively via commands defined
   in DetectorMessenger.
 	
 2- PHYSICS LIST
 
   The PhysicsList contains builders for hadronic interactions.
   Predefined G4 PhysicsConstructors or 'local' PhysicsConstructors can be used 
   (see geant4/source/physics_lists or example runAndEvent/RE04).
   
   In order not to introduce 'artificial' constraints on the step size,
   electromagnetic processes are not registered: there is no continuous energy 
   loss.  
 
   Several hadronic physics options are controlled by environment variables.
   To trigger them, an envHadronic.csh has been added in this example.
   One must select the options wished, and do
        source envHadronic.csh  (or sh) 
 	 
 3- AN EVENT : THE PRIMARY GENERATOR
 
   The primary kinematic consists of a single particle starting at the edge
   of the box. The type of the particle and its energy are set in 
   PrimaryGeneratorAction (neutron 1 MeV), and can be changed via the G4 
   build-in commands of ParticleGun class (see the macros provided with 
   this example).
 	
 4- PHYSICS
 
   An event is killed at the first interaction of the incident paticle.
   The absorption length, also called mean free path, is computed as 
   the mean value of the track length of the incident particle.
   This is why the medium must be 'infinite' : to be sure that interaction
   occurs at any events.
	
   The result is compared with the 'input' value, i.e. with the cross sections
   given by G4HadronicProcessStore and used by Geant4.
   
   The list of nuclear reactions that occured is printed.
   (the number of gamma of deexcitation is not printed).
   
   Then, comes the total list of generated particles and ions.	
   The energy spectrum of the scattered particle (if any) and of the created 
   secondaries are plotted (see SteppingAction).
   
   Momentum conservation is checked as :
   momentum balance = modulus(P_out - P_in)
 	
   A set of macros defining various run conditions are provided.
   The processes can be actived/inactived in order to survey the processes 
   individually.

 5- HISTOGRAMS
         
   The test contains 12 built-in 1D histograms, which are managed by
   G4AnalysisManager and its Messenger. The histos can be individually 
   activated with the command :
   /analysis/h1/set id nbBins  valMin valMax unit 
   where unit is the desired unit for the histo (MeV or keV, etc..)
   (see the macros xxxx.mac).
   
        1	"kinetic energy of scattered primary particle"
	    2	"kinetic energy of gamma"
	    3	"kinetic energy of neutrons"
	    4	"kinetic energy of protons"
	    5	"kinetic energy of deuterons"
	    6	"kinetic energy of alphas"
	    7	"kinetic energy of nuclei"
	    8	"kinetic energy of mesons"
	    9	"kinetic energy of baryons"
	    10	"Q = Ekin out - Ekin in"
	    11	"Pbalance = mag(P_out - P_in)"
	    12	"atomic mass of nuclei"				
      
   The histograms are managed by the HistoManager class and its Messenger. 
   The histos can be individually activated with the command :
   /analysis/h1/set id nbBins  valMin valMax unit 
   where unit is the desired unit for the histo (MeV or keV, deg or mrad, etc..)
   
   One can control the name of the histograms file with the command:
   /analysis/setFileName  name  (default Hadr03)
   
   It is possible to choose the format of the histogram file : root (default),
   xml, csv, by using namespace in HistoManager.hh
       
   It is also possible to print selected histograms on an ascii file:
   /analysis/h1/setAscii id
   All selected histos will be written on a file name.ascii (default Hadr03) 
 	 				
 6- VISUALIZATION
 
   The Visualization Manager is set in the main().
   The initialisation of the drawing is done via the commands
   /vis/... in the macro vis.mac. To get visualisation:
   > /control/execute vis.mac
 	
   The detector has a default view which is a longitudinal view of the box.
   The tracks are drawn at the end of event, and erased at the end of run.
	
 7- HOW TO START ?
 
   Execute Hadr03 in 'batch' mode from macro files :
 	% Hadr03   inelastic.mac
 		
   Execute Hadr03 in 'interactive mode' with visualization :
 	% Hadr03
	Idle> control/execute vis.mac
 	....
 	Idle> type your commands
 	....
 	Idle> exit
