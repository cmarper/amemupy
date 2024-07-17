import ROOT
from sl3offsets import sl3_ofssets
from pprint import pprint

def print_tree(input_file, output_file, station=None):
    f = ROOT.TFile.Open(input_file)
    tree = f.Get("dtNtupleProducer/DTTREE")  # DTNtuples
    with open(output_file, "w") as out:
        for iev, event in enumerate(tree):
        #for iev, event in enumerate(f.dtNtupleProducer.DTTREE):
        #for event in f.DTTree:
            #if iev == 10000: break
            bx = event.event_bunchCrossing  # DTNtuples
            for elem in zip(event.ph2Digi_wheel, event.ph2Digi_sector, event.ph2Digi_station, event.ph2Digi_superLayer, event.ph2Digi_layer, event.ph2Digi_wire, event.ph2Digi_time):
                if station:
                    if elem[2] != station: continue
                shift = sl3_ofssets[(elem[0], elem[1], elem[2])]
                out.write("{} {} {} {} {} {} {} {} {}\n".format(elem[0], elem[1], elem[2], shift / 10., elem[3] - 1, elem[4] - 1, elem[5] - 1, int(elem[6]), bx))
            out.write("-1\n")

def print_emulator(input_file, output_file):
    f = ROOT.TFile.Open(input_file)
    tree = f.Get("dtNtupleProducer/DTTREE")
    with open(output_file, "w") as out:
        for iev, event in enumerate(tree):
            for elem in zip(event.ph2TpgPhiEmuAm_quality, event.ph2TpgPhiEmuAm_posLoc_x, event.ph2TpgPhiEmuAm_dirLoc_phi,
                    event.ph2TpgPhiEmuAm_t0, event.ph2TpgPhiEmuAm_phi, event.ph2TpgPhiEmuAm_phiB, event.ph2TpgPhiEmuAm_chi2,
                    event.ph2TpgPhiEmuAm_wheel, event.ph2TpgPhiEmuAm_sector, event.ph2TpgPhiEmuAm_station,
                    event.ph2TpgPhiEmuAm_superLayer,
                    event.ph2TpgPhiEmuAm_pathWireId, event.ph2TpgPhiEmuAm_pathTDC, event.ph2TpgPhiEmuAm_pathLat):
                out.write(" ".join([str(e) for e in elem[0:11]]))
                for i in range(8):
                    out.write(" %s" % elem[11][i])
                for i in range(8):
                    out.write(" %s" % elem[12][i])
                for i in range(8):
                    out.write(" %s" % elem[13][i])
                out.write("\n")
            out.write("-1\n")


def print_opendata_tree(input_file, output_file, station=None):
    f = ROOT.TFile.Open(input_file)
    tree = f.Get("DTTree")  # OPEN DATA
    with open(output_file, "w") as out:
        for iev, event in enumerate(tree):
        #for iev, event in enumerate(f.dtNtupleProducer.DTTREE):
        #for event in f.DTTree:
            #if iev == 10000: break
            bx = event.bunchXing  # OPEN DATA
            for elem in zip(event.digi_wheel, event.digi_sector, event.digi_station, event.digi_sl, event.digi_layer, event.digi_wire, event.digi_time):
                if station:
                    if elem[2] != station: continue
                shift = sl3_ofssets[(elem[0], elem[1], elem[2])]
                out.write("{} {} {} {} {} {} {} {} {}\n".format(elem[0], elem[1], elem[2], shift / 10., elem[3] - 1, elem[4] - 1, elem[5] - 1, int(elem[6]), bx))
            out.write("-1\n")

# def print_tree_consid_emulator(input_file, output_file, station=None):
    # f = ROOT.TFile.Open(input_file)
    # tree = f.Get("dtNtupleProducer/DTTREE")
    # with open(output_file, "w") as out:
        # for iev, event in enumerate(tree):
        # #for iev, event in enumerate(f.dtNtupleProducer.DTTREE):
        # #for event in f.DTTree:
            # #if iev == 10000: break
            # bx = event.event_bunchCrossing            
            # emu_hits = {}
            # for elem in zip(event.ph2TpgPhiEmuAm_wheel, event.ph2TpgPhiEmuAm_sector, event.ph2TpgPhiEmuAm_station, event.ph2TpgPhiEmuAm_pathWireId, event.ph2TpgPhiEmuAm_pathTDC):
                # if (elem[0], elem[1], elem[2]) not in emu_hits:
                    # emu_hits[(elem[0], elem[1], elem[2])] = []
                
                # for hit in zip(elem[3], elem[4]):
                    # if hit[0] != -1 and hit[1] != -1:
                        # emu_hits[(elem[0], elem[1], elem[2])].append((hit[0], hit[1]))
            
            # pprint(emu_hits)
            
            # for elem in zip(event.ph2Digi_wheel, event.ph2Digi_sector, event.ph2Digi_station, event.ph2Digi_superLayer, event.ph2Digi_layer, event.ph2Digi_wire, event.ph2Digi_time):
                # shift = sl3_ofssets[(elem[0], elem[1], elem[2])]
                # if (elem[0], elem[1], elem[2]) not in emu_hits:
                    # out.write("{} {} {} {} {} {} {} {} {}\n".format(elem[0], elem[1], elem[2], shift / 10., elem[3] - 1, elem[4] - 1, elem[5] - 1, int(elem[6]), bx))
                # else:
                    # printed = False
                    # # if (elem[5] - 1, int(elem[6])) not in emu_hits[(elem[0], elem[1], elem[2])]:
                        # # printed = False
                        # # if (elem[5] - 1, int(round(elem[6])) - 1) in emu_hits[(elem[0], elem[1], elem[2])]:
                            # # printed = True
                            # # print "-1", (elem[5] - 1, int(round(elem[6])) - 1)
                            # # out.write("{} {} {} {} {} {} {} {} {}\n".format(elem[0], elem[1], elem[2], shift / 10., elem[3] - 1, elem[4] - 1, elem[5] - 1, int(round(elem[6])) - 1, bx))
                        # # elif (elem[5] - 1, int(round(elem[6])) + 1) in emu_hits[(elem[0], elem[1], elem[2])]:
                            # # printed = True
                            # # out.write("{} {} {} {} {} {} {} {} {}\n".format(elem[0], elem[1], elem[2], shift / 10., elem[3] - 1, elem[4] - 1, elem[5] - 1, int(round(elem[6])) + 1, bx))
                            # # print "+1", (elem[5] - 1, int(round(elem[6])) + 1)
                    
                    # if not printed:
                        # print "0", (elem[5] - 1, int(elem[6]))
                    # out.write("{} {} {} {} {} {} {} {} {}\n".format(elem[0], elem[1], elem[2], shift / 10., elem[3] - 1, elem[4] - 1, elem[5] - 1, int(elem[6]), bx))
            # out.write("-1\n")
            # return


def debug_tree(input_file, station=None):
    f = ROOT.TFile.Open(input_file)
    #for iev, event in enumerate(f.dtNtupleProducer.DTTREE):
    tree = f.Get("dtNtupleProducer/DTTREE")
    for iev, event in enumerate(tree):
        bx = event.event_bunchCrossing
        #if bx > 3500:
        if iev == 488:
            print(iev, bx)
            for elem in zip(event.ph2TpgPhiEmuAm_wheel, event.ph2TpgPhiEmuAm_sector, event.ph2TpgPhiEmuAm_station, event.ph2TpgPhiEmuAm_quality, event.ph2TpgPhiEmuAm_BX):
                if station:
                    if elem[2] != station: continue
                print(elem)
            print("*****************")
            for elem in zip(event.ph2TpgPhiHw_wheel, event.ph2TpgPhiHw_sector, event.ph2TpgPhiHw_station, event.ph2TpgPhiHw_quality, event.ph2TpgPhiHw_BX):
                if station:
                    if elem[2] != station: continue
                print(elem)
            a = raw_input("continue?")
            if a == "":
                continue

def print_am(input_file):
    f = ROOT.TFile.Open(input_file)
    #for iev, event in enumerate(f.dtNtupleProducer.DTTREE):
    tree = f.Get("dtNtupleProducer/DTTREE")
    for iev, event in enumerate(tree):
        for elem in zip(event.ph2TpgPhiEmuAm_quality, event.ph2TpgPhiEmuAm_posLoc_x_raw, event.ph2TpgPhiEmuAm_dirLoc_phi_raw,
                event.ph2TpgPhiEmuAm_t0, event.ph2TpgPhiEmuAm_phi, event.ph2TpgPhiEmuAm_phiB, event.ph2TpgPhiEmuAm_chi2,
                event.ph2TpgPhiEmuAm_wheel, event.ph2TpgPhiEmuAm_sector, event.ph2TpgPhiEmuAm_station,
                event.ph2TpgPhiEmuAm_superLayer,
                event.ph2TpgPhiEmuAm_pathWireId, event.ph2TpgPhiEmuAm_pathTDC, event.ph2TpgPhiEmuAm_pathLat):
            print(" ".join([str(e) for e in elem[0:11]]), end=" ")
            for i in range(8):
                print(elem[11][i], end=" ") #+ (1 if elem[11][i] != -1 else 0),
            for i in range(8):
                print(elem[12][i], end=" ")
                #print elem[12][i] - (10000 if elem[12][i] != -1 else 0),
            for i in range(8):
                print(elem[13][i], end=" ")
            print()
        print(-1)

def print_am_mixer(input_file):
    f = ROOT.TFile.Open(input_file)
    #for iev, event in enumerate(f.dtNtupleProducer.DTTREE):
    tree = f.Get("dtNtupleProducer/DTTREE")
    for iev, event in enumerate(tree):
        for elem in zip(
                event.ph2TpgPhiEmuAm_wheel, event.ph2TpgPhiEmuAm_sector, event.ph2TpgPhiEmuAm_station,
                event.ph2TpgPhiEmuAm_superLayer,
                event.ph2TpgPhiEmuAm_pathWireId, event.ph2TpgPhiEmuAm_pathTDC):
            print(" ".join([str(e) for e in elem[0:4]]), end=" ")
            for i in range(8):
                print(elem[4][i], end=" ")
            for i in range(8):
                print(elem[5][i], end=" ")
            print()
        print(-1)
            


def debug_trees(input_file1, input_file2, station=None):
    f1 = ROOT.TFile.Open(input_file1)
    f2 = ROOT.TFile.Open(input_file2)
    #for iev, event in enumerate(f.dtNtupleProducer.DTTREE):
    tree1 = f1.Get("dtNtupleProducer/DTTREE")
    tree2 = f2.Get("dtNtupleProducer/DTTREE")
    for iev, (event1, event2) in enumerate(zip(tree1, tree2)):
        for elem in zip(event1.ph2TpgPhiEmuAm_wheel, event1.ph2TpgPhiEmuAm_sector, event1.ph2TpgPhiEmuAm_station, event1.ph2TpgPhiEmuAm_quality, event1.ph2TpgPhiEmuAm_BX):
            if station:
                if elem[2] != station: continue
            print(elem)
        print("*****************")
        for elem in zip(event2.ph2TpgPhiEmuAm_wheel, event2.ph2TpgPhiEmuAm_sector, event2.ph2TpgPhiEmuAm_station, event2.ph2TpgPhiEmuAm_quality, event2.ph2TpgPhiEmuAm_BX):
            if station:
                if elem[2] != station: continue
            print(elem)
        a = raw_input("continue?")
        if a == "":
            continue


def print_segs(input_file, output_file = None, station=None):
    f = ROOT.TFile.Open(input_file)
    #with open(output_file, "w") as out:
    for event in f.dtNtupleProducer.DTTREE:
        #for elem in zip(event.ph2Digi_wheel, event.ph2Digi_sector, event.ph2Digi_station, event.ph2Digi_superLayer, event.ph2Digi_layer, event.ph2Digi_wire, event.ph2Digi_time):
        #    print elem
        #for elem in zip(event.ph2Seg_wheel, event.ph2Seg_sector, event.ph2Seg_station, event.ph2Seg , event.ph2Seg_phiHits_time, event.ph2Seg_phiHits_side):
        for elem in zip(event.ph2Seg_wheel, event.ph2Seg_sector, event.ph2Seg_station,
                event.ph2Seg_phi_t0, event.ph2Seg_dirLoc_x, event.ph2Seg_dirLoc_z,
                event.ph2Seg_posLoc_x_SL1, event.ph2Seg_posLoc_x_SL3, event.ph2Seg_posLoc_x_midPlane,
                event.ph2Seg_phi_nHits):
            if elem[3] <= -999: continue
            if station:
                if elem[2] != station: continue
            print(elem[0], elem[1], elem[2], elem[3], elem[4] / elem[5], elem[6], elem[7], elem[8], elem[9])
        print("-1 " * 8 + "-1") 
        #a = raw_input()
        #if a == "":
        #    continue

def iround(f):
    return int(round(f))

def round4(f):
    return round(f, 4)

def fill_events(input_file):
    from math import acos, cos   
    # matching params
    m_minMuPt = 0
    m_maxMuSegDPhi = 0.2
    m_maxMuSegDEta = 0.3
    m_minSegHits = 4    
    
    events = []
    f = ROOT.TFile.Open(input_file)
    for event in f.dtNtupleProducer.DTTREE:
        ev = {}
        # matches = []
        for elem in zip(event.ph2Digi_wheel, event.ph2Digi_sector, event.ph2Digi_station, event.ph2Digi_superLayer, event.ph2Digi_layer, event.ph2Digi_wire, event.ph2Digi_time):
            key = "({}, {}, {})".format(elem[0], elem[1], elem[2])
            if iround(elem[3]) != 1 and iround(elem[3]) != 3: continue
            if not key in ev:
                ev[key] = {}
                ev[key]["hits"] = []
                ev[key]["segments"] = []
            ev[key]["hits"].append({"sl": iround(elem[3]), "la": iround(elem[4]), "wi": iround(elem[5]), "tdc": iround(elem[6])})

        for elem in zip(event.ph2Seg_wheel, event.ph2Seg_sector, event.ph2Seg_station,
                event.ph2Seg_phi_t0, event.ph2Seg_dirLoc_x, event.ph2Seg_dirLoc_z,
                event.ph2Seg_posLoc_x_SL1, event.ph2Seg_posLoc_x_SL3, event.ph2Seg_posLoc_x_midPlane,
                # segment hits
                event.ph2Seg_phiHits_superLayer, event.ph2Seg_phiHits_layer,
                event.ph2Seg_phiHits_wire, event.ph2Seg_phiHits_time, event.ph2Seg_phiHits_side,
                # matching with gen muon
                event.ph2Seg_posGlb_phi, event.ph2Seg_posGlb_eta,
                ):
            key = "({}, {}, {})".format(elem[0], elem[1], elem[2])
            segment = {"t0": round4(elem[3]), "slope": round4(elem[4] / elem[5]),
                "pos_sl1": round4(elem[6]), "pos_sl3": round4(elem[7]), "pos_midplane": round4(elem[8]),
                "hits": []}
            for hit in zip(elem[9], elem[10], elem[11], elem[12], elem[13]):
                segment["hits"].append({"sl": iround(hit[0]), "la": iround(hit[1]), "wi": iround(hit[2]), "tdc": iround(hit[3]), "lat": iround(hit[4])})
                if {"sl": iround(hit[0]), "la": iround(hit[1]), "wi": iround(hit[2]), "tdc": iround(hit[3])} in ev[key]["hits"]:
                    ev[key]["hits"].remove({"sl": iround(hit[0]), "la": iround(hit[1]), "wi": iround(hit[2]), "tdc": iround(hit[3])})

            # matching with gen muons
            matching_flag = 0
            if len(segment["hits"]) >= m_minSegHits:
                for gen_muon in zip(event.gen_pdgId, event.gen_pt, event.gen_phi, event.gen_eta):
                    if abs(gen_muon[0]) != 13 or gen_muon[1] < m_minMuPt: continue
                    muSegDPhi = abs(acos(cos(gen_muon[2] - elem[14])));
                    muSegDEta = abs(gen_muon[3] - elem[15]);

                    if muSegDPhi < m_maxMuSegDPhi and muSegDEta < m_maxMuSegDEta:
                        # in the efficiency results we only consider the segment with highest number of hits, but here
                        # we will consider matched if the eta and phi reqs are satisfied. To be considered later!
                        matching_flag = 1
                        # sign = elem[0] > 0
                        # if ((-1 + 2 * sign) * elem[1]) not in matches:
                            # matches.append((-1 + 2 * sign) * elem[1])

            segment["gen_match"] = matching_flag
            ev[key]["segments"].append(segment)
        events.append(ev)
        # print matches
    return events

def run(input_file):
    f = ROOT.TFile.Open(input_file)
    for event in f.dtNtupleProducer.DTTREE:
        for elem in zip(event.ph2Seg_wheel, event.ph2Seg_sector, event.ph2Seg_station, event.ph2Seg_phi_vDrift):
            print(elem)
        # a = raw_input()
        # if a == "":
            # continue

if __name__ == "__main__":
    # print_tree("/afs/cern.ch/user/r/redondo/public/opendata/DTNtuple_Run194051.root", "open_hits.txt")
    # events = fill_events("/eos/cms/store/group/dpg_dt/comm_dt/background_33BX_neutron/DTDPGNtuple_11_1_0_patch2_Phase2_Simulation_8muInBarrel_woRPC.root")
    # run("/eos/cms/store/group/dpg_dt/comm_dt/background_33BX_neutron/DTDPGNtuple_11_1_0_patch2_Phase2_Simulation_8muInBarrel_woRPC.root")
    #from pprint import pprint
    #pprint(events)
    # import json
    # with open("events.json", "w") as f:
        # json.dump(events, f, indent = 4)

    #print_segs("/afs/cern.ch/work/j/jleonhol/private/am_reimpl_alv/11_2_pkg/CMSSW_11_2_1_patch2/src/DTDPGAnalysis/DTNtuples/test/DTDPGNtuple_10_6_0_Phase2_Simulation.root")
    #print_segs("/eos/cms/store/group/dpg_dt/comm_dt/commissioning_2021_data/ntuples/ST/DTDPGNtuple_run341582.root", station=4)
    #print_tree("/eos/cms/store/group/dpg_dt/comm_dt/commissioning_2021_data/ntuples/ST/DTDPGNtuple_run341539.root", "hits_run341539_mb4.txt", station=4)
    #debug_tree("/eos/cms/store/group/dpg_dt/comm_dt/commissioning_2021_data/ntuples/ST/DTDPGNtuple_run341539.root", "hits_run341539_mb4.txt", station=4)
    #debug_tree("/eos/user/j/jleonhol/ntuplesST/DTDPGNtuple_run341539.root", "hits_run341539_mb4.txt", station=4)
    #debug_trees("/eos/user/j/jleonhol/simulationSamples/mu_PU200_noRPC_noAgeing_newAnalyzer260521.root",
    #    "/eos/user/j/jleonhol/simulationSamples/mu_PU200_noRPC_noAgeing_20210315_cmssw.root")
    # print_am("/eos/user/j/jleonhol/simulationSamples/rossin_noRPC_noAgeing_ext_newSLFitter_full_12_4_2_v9_simple_noconf_1000ev.root")
    #print_am_mixer("/afs/cern.ch/work/j/jleonhol/private/am_reimpl_alv/new_fitter/CMSSW_12_4_2/src/DTDPGAnalysis/DTNtuples/test/DTDPGNtuple_10_6_0_Phase2_Simulation.root")
    #print_am_mixer("/eos/user/j/jleonhol/simulationSamples/rossin_noRPC_noAgeing_ext_newSLFitter_full_12_4_2_v11_simple_mixer_100ev.root")
    #print_am("/eos/user/j/jleonhol/simulationSamples/rossin_noRPC_noAgeing_ext_newSLFitter_full_12_4_2_v13_simple_1000ev.root")
    #print_segs("/eos/user/j/jleonhol/simulationSamples/rossin_noRPC_noAgeing_ext_newSLFitter_full_12_4_2_v6_noconf_100ev.root")
    #print_tree("/eos/user/j/jleonhol/simulationSamples/rossin_noRPC_noAgeing_ext_newSLFitter_full_12_4_2_v5_noconf.root", "hits_12_4.txt")
    #print_tree("/eos/home-j/jleonhol/simulationSamples/data_from_camilo.root", "hits_data.txt")
    #print_segs("/eos/user/j/jleonhol/simulationSamples/mu_PU200_noRPC_noAgeing_ultraext_noconf.root")
    #print_opendata_tree("/afs/cern.ch/user/r/redondo/public/opendata/DTNtuple_194051_194052.root", "hits_194051_194052.txt")
    #print_opendata_tree("/afs/cern.ch/user/r/redondo/public/opendata/DTTree_Run194115_Files2.root", "hits_194115.txt")
    #print_tree("/eos/home-j/jleonhol/simulationSamples/mu_PU200_noRPC_noAgeing_ultraext.root", "hits_simu_10.txt")
   #debug_tree("/eos/home-j/jleonhol/simulationSamples/mu_PU200_noRPC_noAgeing_ultraext_noconf.root")
    #print_am("/eos/home-j/jleonhol/simulationSamples/mu_PU200_noRPC_noAgeing_ultraext_noconf.root")
    #print_am("/eos/home-j/jleonhol/simulationSamples/data_from_camilo.root")
    #print_am("/afs/cern.ch/user/l/llorente/public/dtntup/ntuple-emulator.root")
    #print_am("/eos/user/l/llorente/TriggerPrimitives/CMSSW_13_3_1/src/DTDPGAnalysis/DTNtuples/test/ntuple-emulator.root")
    print_am("/eos/user/l/llorente/TriggerPrimitives/CMSSW_13_3_1/src/DTDPGAnalysis/DTNtuples/test/prueba.root")
