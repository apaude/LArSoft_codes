//////////////////////////////////////////////////////////////////////////////
// Class:       BeamPDS                                                   ///
// File:        BeamPDS_module.cc                                         ///   
//Description:                                                             /// 
//drift BeamPDS calculation module                                        ///
//dumps all the infomrmation in TTrees which needs to analysed             ///
//contact person:apaudel@phys.ksu.edu                                      ///
//////////////////////////////////////////////////////////////////////////////
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "duneprototypes/Protodune/singlephase/DataUtils/ProtoDUNEDataUtils.h"
#include "lardataobj/RawData/RDTimeStamp.h"

#include "larcore/Geometry/Geometry.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCFlux.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/Simulation/AuxDetSimChannel.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/raw.h"
#include "lardataobj/RawData/BeamInfo.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larcoreobj/SummaryData/POTSummary.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/TrackHitMeta.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/EndPoint2D.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larreco/RecoAlg/TrackMomentumCalculator.h"
#include "lardataobj/AnalysisBase/CosmicTag.h"
#include "lardataobj/AnalysisBase/FlashMatch.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"
#include "lardataobj/RawData/OpDetWaveform.h"

#include "TFile.h"
#include "TTree.h"
#include "TDirectory.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TProfile.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TMath.h"
#include "TGraphErrors.h"
#include "TMinuit.h"
#include "TString.h"
#include "TTimeStamp.h"
#include "TVectorD.h"
#include "TCanvas.h"
#include "TFrame.h"
#include "TLine.h"
#include "TAxis.h"
#include "TTimeStamp.h"

//Root and C++ include
#include "TVector3.h"
#include "TLorentzVector.h"
#include <vector>
#include <fstream>
#include "TPaveStats.h"
#include <iostream>
#include <string>
#include "math.h"
#include "stdio.h"
#include <iterator>
#include "TRandom.h"
#include "TSystem.h"
#include "TInterpreter.h"

//using namespace std;


namespace protoana{
  class BeamPDS : public art::EDAnalyzer {
  public:
    explicit BeamPDS(fhicl::ParameterSet const& pset);
    virtual ~BeamPDS();

    void beginJob();
    void endJob();
    void beginRun(const art::Run& run);
    void analyze(const art::Event& evt);
    void reset();
    
  private:
    ProtoDUNEDataUtils fDataUtils;
    TTree* fEventTree;
    geo::GeometryCore const * fGeometry;

    //These are the tree variables I will be using
    Int_t    run;                  
    Int_t    subrun;               
    Int_t    event;
    Short_t Trigger_ID;  
    Double_t evttime; 
    int fNactivefembs[6];
    std::vector<float> trackthetaxz;
    std::vector<float>  trackthetayz;
    std::vector<float> trkstartx;
    std::vector<float> trkstarty;
    std::vector<float> trkstartz;
    std::vector<float> trkendx;
    std::vector<float> trkendy;
    std::vector<float> trkendz;

    std::vector<float> trkstartx_crt2;
    std::vector<float> trkendx_crt2;
    std::vector<float> crtreco_x0;
    std::vector<float> crtreco_x1;
    std::vector<float> crtreco_y0;
    std::vector<float> crtreco_y1;
    std::vector<float> crtreco_z0;
    std::vector<float> crtreco_z1;
    ////photon detector information
    std::vector<short> DAQch_ID;            
    std::vector< std::vector<float> > onda_ID;

    //photon det info
    std::vector<float> trklen;
    std::vector<int> TrkID; 
    std::vector<float>  xprojectedlen;
    //  std::vector<double> t0crt1;
    std::vector<double> t0crt2;
    std::vector<double> crt2tickoffset;
    std::vector<int> tot_trks;
   
    std::string fHitsModuleLabel;
    std::string fTrackModuleLabel;
    std::string fCalorimetryModuleLabel;
    bool  fSaveTrackInfo;
    bool  fSaveCaloInfo;
  };

  //========================================================================
  BeamPDS::BeamPDS(fhicl::ParameterSet const& pset) :
    EDAnalyzer(pset),
    fDataUtils                  (pset.get<fhicl::ParameterSet>("DataUtils")),
    fHitsModuleLabel          (pset.get< std::string >("HitsModuleLabel","")         ), 
    fTrackModuleLabel         (pset.get< std::string >("TrackModuleLabel","")        ),
    fCalorimetryModuleLabel   (pset.get< std::string >("CalorimetryModuleLabel","")  ),
    fSaveTrackInfo            (pset.get< bool>("SaveTrackInfo",false)),
    fSaveCaloInfo             (pset.get< bool>("SaveCaloInfo",false))

  {
    if (fSaveTrackInfo == false) fSaveCaloInfo = false;
  }
 
  //========================================================================
  BeamPDS::~BeamPDS(){
  }
  //========================================================================

  //========================================================================
  void BeamPDS::beginJob(){
    std::cout<<"job begin..."<<std::endl;
    art::ServiceHandle<art::TFileService> tfs;
    fEventTree = tfs->make<TTree>("Event", "Event Tree from Reco");
    fEventTree->Branch("event", &event,"event/I");
    fEventTree->Branch("evttime",&evttime,"evttime/D");
    fEventTree->Branch("run", &run,"run/I");
    fEventTree->Branch("Nactivefembs",&fNactivefembs,"Nactivefembs[6]/I");
    fEventTree->Branch("subrun", &subrun,"surbrun/I");
    fEventTree->Branch("Trigger_ID", &Trigger_ID,"Trigger_ID/S");
    fEventTree->Branch("xprojectedlen",&xprojectedlen);
    fEventTree->Branch("trackthetaxz",&trackthetaxz);
    fEventTree->Branch("trackthetayz",&trackthetayz);
    fEventTree->Branch("trkstartx",&trkstartx);
    fEventTree->Branch("trkstarty",&trkstarty);
    fEventTree->Branch("trkstartz",&trkstartz);
    fEventTree->Branch("trkendx",&trkendx);
    fEventTree->Branch("trkendy",&trkendy);
    fEventTree->Branch("trkendz",&trkendz);
    fEventTree->Branch("trkendx_crt2",&trkendx_crt2);
    fEventTree->Branch("trkstartx_crt2",&trkstartx_crt2);
    fEventTree->Branch("crtreco_x0",&crtreco_x0);
    fEventTree->Branch("crtreco_x1",&crtreco_x1);
    fEventTree->Branch("crtreco_y0",&crtreco_y0);
    fEventTree->Branch("crtreco_y1",&crtreco_y1);
    fEventTree->Branch("crtreco_z0",&crtreco_z0);
    fEventTree->Branch("crtreco_z1",&crtreco_z1);
    fEventTree->Branch("crt2tickoffset",&crt2tickoffset);
    fEventTree->Branch("trklen",&trklen);
    fEventTree->Branch("TrkID",&TrkID);
    fEventTree->Branch("tot_trks",&tot_trks);
    fEventTree->Branch("t0crt2",&t0crt2);
    fEventTree->Branch("DAQch_ID",&DAQch_ID);
    fEventTree->Branch("onda_ID",&onda_ID);
  }

  //========================================================================
  void BeamPDS::endJob(){     

  }

  //========================================================================
  void BeamPDS::beginRun(const art::Run&){
    mf::LogInfo("BeamPDS")<<"begin run..."<<std::endl;
  }
  //========================================================================

  //========================================================================

  //========================================================================

  void BeamPDS::analyze( const art::Event& evt){//analyze
    reset();  

    // art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
    fGeometry = &*(art::ServiceHandle<geo::Geometry>());
    //Detector properties service
    auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService>()->DataFor(evt);
   
    art::Handle< std::vector<recob::Track> > trackListHandle;
    std::vector<art::Ptr<recob::Track> > tracklist;
    if(evt.getByLabel("pandoraTrack",trackListHandle)){
      art::fill_ptr_vector(tracklist, trackListHandle);
    }
    else return;


    art::Handle< std::vector<recob::PFParticle> > PFPListHandle; 
    std::vector<art::Ptr<recob::PFParticle> > pfplist;

    if(evt.getByLabel("pandora",PFPListHandle)) art::fill_ptr_vector(pfplist, PFPListHandle);
  
   
    art::Handle< std::vector<recob::Hit> > hitListHandle; // to get information about the hits
    std::vector<art::Ptr<recob::Hit>> hitlist;
    if(evt.getByLabel(fHitsModuleLabel, hitListHandle))
      art::fill_ptr_vector(hitlist, hitListHandle);
    art::FindManyP<recob::Hit, recob::TrackHitMeta> fmthm(trackListHandle, evt, fTrackModuleLabel); // to associate tracks and hits
    art::FindManyP<anab::T0> trk_t0_assn_v(PFPListHandle, evt ,"pandora");
    art::FindManyP<recob::PFParticle> pfp_trk_assn(trackListHandle,evt,"pandoraTrack");
    art::FindManyP<anab::T0> fmT0(trackListHandle, evt ,"pmtrack");
    std::vector<const sim::SimChannel*> fSimChannels;
    try{
      evt.getView("largeant", fSimChannels);
    }catch (art::Exception const&e){
    }

    //Get 2-CRT T0
    art::FindManyP<anab::T0> fmt0crt2(trackListHandle, evt, "crtreco");
    art::FindManyP<anab::CosmicTag> fmctcrt2(trackListHandle, evt, "crtreco");
  
    //photon trigger info
    auto const TriggHandle = evt.getValidHandle<std::vector<raw::RDTimeStamp>>("timingrawdecoder:daq");
    std::vector<raw::RDTimeStamp> TriggVec(*TriggHandle);
    raw::RDTimeStamp theTrigger = TriggVec[0];
    const double & RDTSTrigger = theTrigger.GetFlags();

    /////////////////



    run = evt.run();
    Trigger_ID=RDTSTrigger;
    subrun = evt.subRun();
    event = evt.id().event();
    art::Timestamp ts = evt.time();
    TTimeStamp tts(ts.timeHigh(), ts.timeLow());
    evttime=tts.AsDouble();



    // Get number of active fembs
    if(!evt.isRealData()){
      for(int k=0; k < 6; k++)
	fNactivefembs[k] = 20;
    }
    else{
      for(int k=0; k < 6; k++)
	fNactivefembs[k] = fDataUtils.GetNActiveFembsForAPA(evt, k);
    }


    int trks=0;
    size_t NTracks = tracklist.size();
    int entries_count=0;
    for(size_t i=0; i<NTracks;++i){
      double xoffset=0.0;
      art::Ptr<recob::Track> ptrack(trackListHandle, i);
      const recob::Track& track = *ptrack;

      double this_t0crt2=-DBL_MAX;
      //double this_t0crt1=-DBL_MAX;
      double this_crt2x0 = -DBL_MAX;
      double this_crt2x1 = -DBL_MAX;
      double this_crt2y0 = -DBL_MAX;
      double this_crt2y1 = -DBL_MAX;
      double this_crt2z0 = -DBL_MAX;
      double this_crt2z1 = -DBL_MAX;

      bool test1=true;
      if(fmt0crt2.isValid()){
	auto const& vt0crt2 = fmt0crt2.at(i);
	if (!vt0crt2.empty()){
	  this_t0crt2 = vt0crt2[0]->Time();
	  test1=false;
	}
	
      }
      if(test1) continue;
      if (fmctcrt2.isValid()){
	auto const& vctcrt2 = fmctcrt2.at(i);
	if (!vctcrt2.empty()){
	  this_crt2x0 = vctcrt2[0]->EndPoint1()[0];
	  this_crt2x1 = vctcrt2[0]->EndPoint2()[0];
	  this_crt2y0 = vctcrt2[0]->EndPoint1()[1];
	  this_crt2y1 = vctcrt2[0]->EndPoint2()[1];
	  this_crt2z0 = vctcrt2[0]->EndPoint1()[2];
	  this_crt2z1 = vctcrt2[0]->EndPoint2()[2];

      }
    }
      

      ///////////Storing photon detector signals/////////////////////
      entries_count++;
      short aai=0; short DAQch=0;
      if(entries_count==1){
	auto const& waveforms = *evt.getValidHandle<std::vector<raw::OpDetWaveform>>("ssprawdecoder:external");

	for (auto const& in_wave : waveforms){
	  DAQch=in_wave.ChannelNumber();
	  DAQch_ID.push_back(DAQch); 
	  onda_ID.push_back(std::vector<float>());
	  for(short i=0; i<2000; i++){
	    aai=in_wave[i]; 
	    onda_ID.back().push_back(aai);
	  } 
	}
      }
      //////////////////////////////////////////////////////////////
      auto pos = track.Vertex();
      auto end = track.End();
      auto dir_start = track.VertexDirection();

      double theta_xz = std::atan2(dir_start.X(), dir_start.Z());
      double theta_yz = std::atan2(dir_start.Y(), dir_start.Z());
      auto allHits=fmthm.at(i);
      double ticksoffset=0;
      if (this_t0crt2 > -DBL_MAX) ticksoffset = this_t0crt2/500.+detProp.GetXTicksOffset (allHits[0]->WireID().Plane, allHits[0]->WireID().TPC, allHits[0]->WireID().Cryostat);
      xoffset = detProp.ConvertTicksToX(ticksoffset,allHits[0]->WireID());
      std::cout<<"tickoffset , x offset "<<ticksoffset<<"  "<<xoffset<<" default term "<<detProp.GetXTicksOffset (allHits[10]->WireID().Plane, allHits[10]->WireID().TPC, allHits[10]->WireID().Cryostat)<<std::endl;

      trks++;
      xprojectedlen.push_back(TMath::Abs(end.X()-pos.X()));
      trackthetaxz.push_back(theta_xz);
      trackthetayz.push_back(theta_yz);
      trkstartx.push_back(pos.X());
      trkstartx_crt2.push_back(pos.X()-xoffset);
      trkstarty.push_back(pos.Y());
      trkstartz.push_back(pos.Z());
      trkendx.push_back(end.X());
      trkendx_crt2.push_back(end.X()-xoffset);
      crtreco_x0.push_back(this_crt2x0);
      crtreco_x1.push_back(this_crt2x1);
      crtreco_y0.push_back(this_crt2y0);
      crtreco_y1.push_back(this_crt2y1);
      crtreco_z0.push_back(this_crt2z0);
      crtreco_z1.push_back(this_crt2z1);

      trkendy.push_back(end.Y());
      trkendz.push_back(end.Z());
      trklen.push_back(track.Length());
      TrkID.push_back(track.ID());
      crt2tickoffset.push_back(ticksoffset);
      t0crt2.push_back(this_t0crt2/500.0);
    } // loop over trks...
    tot_trks.push_back(trks);
    fEventTree->Fill();
  } // end of analyze function
	   
  /////////////////// Defintion of reset function ///////////
  void BeamPDS::reset(){
    run = -9999;
    subrun = -9999;
    event = -9999;
    evttime = -9999;
    //all_trks = -9999;
    for(int k=0; k < 6; k++)
      fNactivefembs[k] = -9999;
    trackthetaxz.clear();
    trackthetayz.clear();
    trkstartx.clear();
    trkstartx_crt2.clear();
    trkendx_crt2.clear();
    trkstarty.clear();
    trkstartz.clear();
  
    trkendx.clear();
    trkendy.clear();
    trkendz.clear();
    trklen.clear();
    TrkID.clear();
    tot_trks.clear();
    xprojectedlen.clear();
    crtreco_x0.clear();
    crtreco_x1.clear();
    crtreco_y0.clear();
    crtreco_y1.clear();
    crtreco_z0.clear();
    crtreco_z1.clear();
    crt2tickoffset.clear();
    DAQch_ID.clear();
    onda_ID.clear();
    //  t0crt1.clear();
    t0crt2.clear();

  }
  //////////////////////// End of definition ///////////////	
	  
  DEFINE_ART_MODULE(BeamPDS)
}
