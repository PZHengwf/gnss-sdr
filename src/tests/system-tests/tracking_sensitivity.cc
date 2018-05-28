/*!
 * \file tracking_sensitivity.cc
 * \brief  This class implements a test for measuring tracking sensitivity
 * \author Carles Fernandez-Prades, 2018. cfernandez(at)cttc.es
 *
 *
 * -------------------------------------------------------------------------
 *
 * Copyright (C) 2010-2018  (see AUTHORS file for a list of contributors)
 *
 * GNSS-SDR is a software defined Global Navigation
 *          Satellite Systems receiver
 *
 * This file is part of GNSS-SDR.
 *
 * GNSS-SDR is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * GNSS-SDR is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with GNSS-SDR. If not, see <http://www.gnu.org/licenses/>.
 *
 * -------------------------------------------------------------------------
 */

#include "concurrent_map.h"
#include "concurrent_queue.h"
#include "control_thread.h"
#include "in_memory_configuration.h"
#include "file_configuration.h"
#include "MATH_CONSTANTS.h"
#include "gnuplot_i.h"
#include "test_flags.h"
#include "signal_generator_flags.h"
#include "tracking_dump_reader.h"
#include "armadillo"
#include <boost/filesystem.hpp>
#include <glog/logging.h>
#include <gpstk/RinexUtilities.hpp>
#include <gpstk/Rinex3ObsBase.hpp>
#include <gpstk/Rinex3ObsData.hpp>
#include <gpstk/Rinex3ObsHeader.hpp>
#include <gpstk/Rinex3ObsStream.hpp>
#include <gtest/gtest.h>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <numeric>
#include <thread>

DEFINE_string(config_file_trk, std::string(""), "File containing the configuration parameters for the Tracking Sensitivity test.");
DEFINE_int32(fs_in, 4000000, "Sampling rate, in Samples/s");

concurrent_queue<Gps_Acq_Assist> global_gps_acq_assist_queue;
concurrent_map<Gps_Acq_Assist> global_gps_acq_assist_map;

class TrackingSensitivityTest : public ::testing::Test
{
public:
    void config_1();
    void config_2();
    bool check_valid_rinex_obs(std::string filename, int rinex_ver);
    void read_rinex_files(
        std::vector<arma::mat>& pseudorange_meas,
        std::vector<arma::mat>& carrierphase_meas,
        std::vector<arma::mat>& doppler_meas,
        std::vector<arma::mat>& cn0_meas,
        arma::mat& sow_prn_ref,
        int signal_type);
    void read_dump_files();

    std::shared_ptr<InMemoryConfiguration> config;
    std::shared_ptr<FileConfiguration> config2;

    std::string filename_rinex_obs = FLAGS_config_file_trk;  ////

    const double central_freq = 1575420000.0;
    const double gain_dB = 40.0;

    const int number_of_taps = 11;
    const int number_of_bands = 2;
    const float band1_begin = 0.0;
    const float band1_end = 0.48;
    const float band2_begin = 0.52;
    const float band2_end = 1.0;
    const float ampl1_begin = 1.0;
    const float ampl1_end = 1.0;
    const float ampl2_begin = 0.0;
    const float ampl2_end = 0.0;
    const float band1_error = 1.0;
    const float band2_error = 1.0;
    const int grid_density = 16;
    const float zero = 0.0;
    const int number_of_channels = 8;
    const int in_acquisition = 1;

    const float threshold = 0.01;
    const float doppler_max = 8000.0;
    const float doppler_step = 500.0;
    const int max_dwells = 1;
    const int tong_init_val = 2;
    const int tong_max_val = 10;
    const int tong_max_dwells = 30;
    const int coherent_integration_time_ms = 1;

    const float pll_bw_hz = 30.0;
    const float dll_bw_hz = 4.0;
    const float early_late_space_chips = 0.5;

    const int display_rate_ms = 500;
    const int output_rate_ms = 100;
    const int averaging_depth = 10;

    const int num_prn_gps = 33;
    const int num_prn_gal = 31;
};


void TrackingSensitivityTest::config_1()
{
    config = std::make_shared<InMemoryConfiguration>();

    config->set_property("GNSS-SDR.internal_fs_sps", std::to_string(FLAGS_fs_in));

    // Set the assistance system parameters
    config->set_property("GNSS-SDR.SUPL_gps_ephemeris_server", "supl.google.com");
    config->set_property("GNSS-SDR.SUPL_gps_ephemeris_port", std::to_string(7275));
    config->set_property("GNSS-SDR.SUPL_gps_acquisition_server", "supl.google.com");
    config->set_property("GNSS-SDR.SUPL_gps_acquisition_port", std::to_string(7275));
    config->set_property("GNSS-SDR.SUPL_MCC", std::to_string(244));
    config->set_property("GNSS-SDR.SUPL_MNS", std::to_string(5));
    config->set_property("GNSS-SDR.SUPL_LAC", "0x59e2");
    config->set_property("GNSS-SDR.SUPL_CI", "0x31b0");

    // Set the Signal Source
    config->set_property("SignalSource.item_type", "cshort");
    config->set_property("SignalSource.implementation", "UHD_Signal_Source");
    config->set_property("SignalSource.freq", std::to_string(central_freq));
    config->set_property("SignalSource.sampling_frequency", std::to_string(FLAGS_fs_in));
    config->set_property("SignalSource.gain", std::to_string(gain_dB));
    //config->set_property("SignalSource.subdevice", FLAGS_subdevice);
    config->set_property("SignalSource.samples", std::to_string(FLAGS_fs_in * FLAGS_duration));
    //config->set_property("SignalSource.device_address", FLAGS_device_address);

    // Set the Signal Conditioner
    config->set_property("SignalConditioner.implementation", "Signal_Conditioner");
    config->set_property("DataTypeAdapter.implementation", "Pass_Through");
    config->set_property("DataTypeAdapter.item_type", "cshort");
    config->set_property("InputFilter.implementation", "Fir_Filter");
    config->set_property("InputFilter.dump", "false");
    config->set_property("InputFilter.input_item_type", "cshort");
    config->set_property("InputFilter.output_item_type", "gr_complex");
    config->set_property("InputFilter.taps_item_type", "float");
    config->set_property("InputFilter.number_of_taps", std::to_string(number_of_taps));
    config->set_property("InputFilter.number_of_bands", std::to_string(number_of_bands));
    config->set_property("InputFilter.band1_begin", std::to_string(band1_begin));
    config->set_property("InputFilter.band1_end", std::to_string(band1_end));
    config->set_property("InputFilter.band2_begin", std::to_string(band2_begin));
    config->set_property("InputFilter.band2_end", std::to_string(band2_end));
    config->set_property("InputFilter.ampl1_begin", std::to_string(ampl1_begin));
    config->set_property("InputFilter.ampl1_end", std::to_string(ampl1_end));
    config->set_property("InputFilter.ampl2_begin", std::to_string(ampl2_begin));
    config->set_property("InputFilter.ampl2_end", std::to_string(ampl2_end));
    config->set_property("InputFilter.band1_error", std::to_string(band1_error));
    config->set_property("InputFilter.band2_error", std::to_string(band2_error));
    config->set_property("InputFilter.filter_type", "bandpass");
    config->set_property("InputFilter.grid_density", std::to_string(grid_density));
    config->set_property("InputFilter.sampling_frequency", std::to_string(FLAGS_fs_in));
    config->set_property("InputFilter.IF", std::to_string(zero));
    config->set_property("Resampler.implementation", "Pass_Through");
    config->set_property("Resampler.dump", "false");
    config->set_property("Resampler.item_type", "gr_complex");
    config->set_property("Resampler.sample_freq_in", std::to_string(FLAGS_fs_in));
    config->set_property("Resampler.sample_freq_out", std::to_string(FLAGS_fs_in));

    // Set the number of Channels
    config->set_property("Channels_1C.count", std::to_string(number_of_channels));
    config->set_property("Channels.in_acquisition", std::to_string(in_acquisition));
    config->set_property("Channel.signal", "1C");

    // Set Acquisition
    config->set_property("Acquisition_1C.implementation", "GPS_L1_CA_PCPS_Tong_Acquisition");
    config->set_property("Acquisition_1C.item_type", "gr_complex");
    config->set_property("Acquisition_1C.if", std::to_string(zero));
    config->set_property("Acquisition_1C.coherent_integration_time_ms", std::to_string(coherent_integration_time_ms));
    config->set_property("Acquisition_1C.threshold", std::to_string(threshold));
    config->set_property("Acquisition_1C.doppler_max", std::to_string(doppler_max));
    config->set_property("Acquisition_1C.doppler_step", std::to_string(doppler_step));
    config->set_property("Acquisition_1C.bit_transition_flag", "false");
    config->set_property("Acquisition_1C.max_dwells", std::to_string(max_dwells));
    config->set_property("Acquisition_1C.tong_init_val", std::to_string(tong_init_val));
    config->set_property("Acquisition_1C.tong_max_val", std::to_string(tong_max_val));
    config->set_property("Acquisition_1C.tong_max_dwells", std::to_string(tong_max_dwells));

    // Set Tracking
    config->set_property("Tracking_1C.implementation", "GPS_L1_CA_DLL_PLL_Tracking");
    config->set_property("Tracking_1C.item_type", "gr_complex");
    config->set_property("Tracking_1C.if", std::to_string(zero));
    config->set_property("Tracking_1C.dump", "false");
    config->set_property("Tracking_1C.dump_filename", "./tracking_ch_");
    config->set_property("Tracking_1C.pll_bw_hz", std::to_string(pll_bw_hz));
    config->set_property("Tracking_1C.dll_bw_hz", std::to_string(dll_bw_hz));
    config->set_property("Tracking_1C.early_late_space_chips", std::to_string(early_late_space_chips));

    // Set Telemetry
    config->set_property("TelemetryDecoder_1C.implementation", "GPS_L1_CA_Telemetry_Decoder");
    config->set_property("TelemetryDecoder_1C.dump", "false");

    // Set Observables
    config->set_property("Observables.implementation", "Hybrid_Observables");
    config->set_property("Observables.dump", "false");
    config->set_property("Observables.dump_filename", "./observables.dat");

    // Set PVT
    config->set_property("PVT.implementation", "RTKLIB_PVT");
    config->set_property("PVT.output_rate_ms", std::to_string(output_rate_ms));
    config->set_property("PVT.display_rate_ms", std::to_string(display_rate_ms));
    config->set_property("PVT.dump_filename", "./PVT");
    config->set_property("PVT.nmea_dump_filename", "./gnss_sdr_pvt.nmea");
    config->set_property("PVT.flag_nmea_tty_port", "false");
    config->set_property("PVT.nmea_dump_devname", "/dev/pts/4");
    config->set_property("PVT.flag_rtcm_server", "false");
    config->set_property("PVT.flag_rtcm_tty_port", "false");
    config->set_property("PVT.rtcm_dump_devname", "/dev/pts/1");
    config->set_property("PVT.dump", "false");
}


void TrackingSensitivityTest::config_2()
{
    if (FLAGS_config_file_trk.empty())
        {
            std::string path = std::string(TEST_PATH);
            std::string filename = path + "../../conf/gnss-sdr_GPS_L1_USRP_X300_realtime.conf";
            config2 = std::make_shared<FileConfiguration>(filename);
        }
    else
        {
            config2 = std::make_shared<FileConfiguration>(FLAGS_config_file_trk);
        }

    int d_sampling_rate;
    d_sampling_rate = config2->property("GNSS-SDR.internal_fs_sps", FLAGS_fs_in);
    config2->set_property("SignalSource.samples", std::to_string(d_sampling_rate * FLAGS_duration));
}


bool TrackingSensitivityTest::check_valid_rinex_obs(std::string filename, int rinex_ver)
{
    bool res = false;
    if (rinex_ver == 2)
        {
            res = gpstk::isRinexObsFile(filename);
        }
    if (rinex_ver == 3)
        {
            res = gpstk::isRinex3ObsFile(filename);
        }
    return res;
}


void TrackingSensitivityTest::read_rinex_files(
    std::vector<arma::mat>& pseudorange_meas,
    std::vector<arma::mat>& carrierphase_meas,
    std::vector<arma::mat>& doppler_meas,
    std::vector<arma::mat>& cn0_meas,
    arma::mat& sow_prn_ref,
    int signal_type)
{
    bool ref_exist = false;
    bool meas_exist = false;
    gpstk::SatID::SatelliteSystem sat_type = gpstk::SatID::systemUnknown;
    int max_prn = 0;
    std::string pr_string;
    std::string cp_string;
    std::string dp_string;
    std::string cn0_string;
    std::string signal_type_string;
    sow_prn_ref.reset();

    switch (signal_type)
        {
        case 0:  //GPS L1

            sat_type = gpstk::SatID::systemGPS;
            max_prn = num_prn_gps;
            pr_string = "C1C";
            cp_string = "L1C";
            dp_string = "D1C";
            cn0_string = "S1C";
            signal_type_string = "GPS L1 C/A";
            break;

        case 1:  //Galileo E1B

            sat_type = gpstk::SatID::systemGalileo;
            max_prn = num_prn_gal;
            pr_string = "C1B";
            cp_string = "L1B";
            dp_string = "D1B";
            cn0_string = "S1B";
            signal_type_string = "Galileo E1B";
            break;

        case 2:  //GPS L5

            sat_type = gpstk::SatID::systemGPS;
            max_prn = num_prn_gps;
            pr_string = "C5X";
            cp_string = "L5X";
            dp_string = "D5X";
            cn0_string = "S5X";
            signal_type_string = "GPS L5";
            break;

        case 3:  //Galileo E5a

            sat_type = gpstk::SatID::systemGalileo;
            max_prn = num_prn_gal;
            pr_string = "C5X";
            cp_string = "L5X";
            dp_string = "D5X";
            cn0_string = "S5X";
            signal_type_string = "Galileo E5a";
            break;
        }

    // Open and read reference RINEX observables file
    std::cout << "Read: RINEX " << signal_type_string << " True" << std::endl;
    try
        {
            gpstk::Rinex3ObsStream r_ref(filename_rinex_obs);
            r_ref.exceptions(std::ios::failbit);
            gpstk::Rinex3ObsData r_ref_data;
            gpstk::Rinex3ObsHeader r_ref_header;
            gpstk::RinexDatum dataobj;
            r_ref >> r_ref_header;

            while (r_ref >> r_ref_data)
                {
                    for (int myprn = 1; myprn < max_prn; myprn++)
                        {
                            gpstk::SatID prn(myprn, sat_type);
                            gpstk::CommonTime time = r_ref_data.time;
                            double sow(static_cast<gpstk::GPSWeekSecond>(time).sow);
                            gpstk::Rinex3ObsData::DataMap::iterator pointer = r_ref_data.obs.find(prn);
                            if (pointer == r_ref_data.obs.end())
                                {
                                    // PRN not present; do nothing
                                }
                            else
                                {
                                    dataobj = r_ref_data.getObs(prn, pr_string, r_ref_header);
                                    double P1 = dataobj.data;
                                    pseudorange_meas.at(myprn).insert_rows(pseudorange_meas.at(myprn).n_rows, arma::rowvec({sow, P1}));

                                    dataobj = r_ref_data.getObs(prn, cp_string, r_ref_header);
                                    double L1 = dataobj.data;
                                    carrierphase_meas.at(myprn).insert_rows(carrierphase_meas.at(myprn).n_rows, arma::rowvec({sow, L1}));

                                    dataobj = r_ref_data.getObs(prn, dp_string, r_ref_header);
                                    double D1 = dataobj.data;
                                    doppler_meas.at(myprn).insert_rows(doppler_meas.at(myprn).n_rows, arma::rowvec({sow, D1}));

                                    dataobj = r_ref_data.getObs(prn, cn0_string, r_ref_header);
                                    double S1 = dataobj.data;
                                    cn0_meas.at(myprn).insert_rows(cn0_meas.at(myprn).n_rows, arma::rowvec({sow, S1}));

                                    ref_exist = true;
                                }  // End of 'if( pointer == roe.obs.end() )'
                        }          // end for
                }                  // end while
        }                          // End of 'try' block
    catch (const gpstk::FFStreamError& e)
        {
            std::cout << e;
            exit(1);
        }
    catch (const gpstk::Exception& e)
        {
            std::cout << e;
            exit(1);
        }
    catch (...)
        {
            std::cout << "unknown error.  I don't feel so well..." << std::endl;
            exit(1);
        }
}


void TrackingSensitivityTest::read_dump_files()
{
    //load the measured values
    tracking_dump_reader trk_dump;

    ASSERT_EQ(trk_dump.open_obs_file(std::string("./tracking_ch_0.dat")), true)
        << "Failure opening tracking dump file";

    long int nepoch = trk_dump.num_epochs();
    std::cout << "Measured observation epochs=" << nepoch << std::endl;

    arma::vec trk_timestamp_s = arma::zeros(nepoch, 1);
    //arma::vec trk_acc_carrier_phase_cycles = arma::zeros(nepoch, 1);
    //arma::vec trk_Doppler_Hz = arma::zeros(nepoch, 1);
    //arma::vec trk_prn_delay_chips = arma::zeros(nepoch, 1);

    //std::vector<double> prompt;
    //std::vector<double> early;
    //std::vector<double> late;
    //std::vector<double> promptI;
    //std::vector<double> promptQ;
    std::vector<double> CN0;
    std::vector<unsigned int> PRN;

    int epoch_counter = 0;
    while (trk_dump.read_binary_obs())
        {
            trk_timestamp_s(epoch_counter) = static_cast<double>(trk_dump.PRN_start_sample_count) / static_cast<double>(FLAGS_fs_in);
            //trk_acc_carrier_phase_cycles(epoch_counter) = trk_dump.acc_carrier_phase_rad / GPS_TWO_PI;
            //trk_Doppler_Hz(epoch_counter) = trk_dump.carrier_doppler_hz;

            //double delay_chips = GPS_L1_CA_CODE_LENGTH_CHIPS - GPS_L1_CA_CODE_LENGTH_CHIPS * (fmod((static_cast<double>(trk_dump.PRN_start_sample_count) + trk_dump.aux1) / static_cast<double>(baseband_sampling_freq), 1.0e-3) / 1.0e-3);

            //trk_prn_delay_chips(epoch_counter) = delay_chips;
            epoch_counter++;
            //prompt.push_back(trk_dump.abs_P);
            //early.push_back(trk_dump.abs_E);
            //late.push_back(trk_dump.abs_L);
            //promptI.push_back(trk_dump.prompt_I);
            //promptQ.push_back(trk_dump.prompt_Q);
            CN0.push_back(trk_dump.CN0_SNV_dB_Hz);
            PRN.push_back(trk_dump.PRN);
        }
}


TEST_F(TrackingSensitivityTest, GPSL1)
{
    // Create a new ControlThread object with a smart pointer
    std::shared_ptr<ControlThread> control_thread;
    if (FLAGS_config_file_trk.empty())
        {
            control_thread = std::make_shared<ControlThread>(config);
        }
    else
        {
            control_thread = std::make_shared<ControlThread>(config2);
        }

    // Run receiver
    try
        {
            control_thread->run();
        }
    catch (const boost::exception& e)
        {
            std::cout << "Boost exception: " << boost::diagnostic_information(e);
        }
    catch (const std::exception& ex)
        {
            std::cout << "STD exception: " << ex.what();
        }

    // Read generated RINEX file
    std::cout << "Validating generated RINEX obs file: " << filename_rinex_obs << " ..." << std::endl;
    bool is_rinex_obs_valid = check_valid_rinex_obs(filename_rinex_obs, 3);
    ASSERT_EQ(true, is_rinex_obs_valid) << "The RINEX observation file " << filename_rinex_obs << " is not well formed. Only RINEX v. 3.00 files are allowed";
    std::cout << "The file is valid." << std::endl;

    std::vector<arma::mat> pseudorange_meas(num_prn_gps);
    std::vector<arma::mat> carrierphase_meas(num_prn_gps);
    std::vector<arma::mat> doppler_meas(num_prn_gps);
    std::vector<arma::mat> cn0_meas(num_prn_gps);
    arma::mat sow_prn_ref;
    read_rinex_files(pseudorange_meas, carrierphase_meas, doppler_meas, cn0_meas, sow_prn_ref, 0);
    // Compute minimum tracked CN0
}


int main(int argc, char** argv)
{
    std::cout << "Running Tracking Sensitivity test..." << std::endl;
    int res = 0;
    try
        {
            testing::InitGoogleTest(&argc, argv);
        }
    catch (...)
        {
        }  // catch the "testing::internal::<unnamed>::ClassUniqueToAlwaysTrue" from gtest

    google::ParseCommandLineFlags(&argc, &argv, true);
    google::InitGoogleLogging(argv[0]);

    // Run the Tests
    try
        {
            res = RUN_ALL_TESTS();
        }
    catch (...)
        {
            LOG(WARNING) << "Unexpected catch";
        }

    google::ShutDownCommandLineFlags();
    return res;
}
