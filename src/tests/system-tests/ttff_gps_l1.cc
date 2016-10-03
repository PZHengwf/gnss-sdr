/*!
 * \file ttff_gps_l1.cc
 * \brief  This class implements a test for measuring
 * the Time-To-First-Fix
 * \author Carles Fernandez-Prades, 2016. cfernandez(at)cttc.es
 *
 *
 * -------------------------------------------------------------------------
 *
 * Copyright (C) 2010-2016  (see AUTHORS file for a list of contributors)
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

#include <chrono>
#include <cstdlib>
#include <ctime>
#include <cerrno>
#include <numeric>
#include <string>
#include <sys/types.h>
#include <sys/ipc.h>
#include <sys/msg.h>
#include <thread>
#include <gflags/gflags.h>
#include <glog/logging.h>
#include <gtest/gtest.h>
#include "concurrent_map.h"
#include "concurrent_queue.h"
#include "control_thread.h"
#include "file_configuration.h"
#include "in_memory_configuration.h"
#include "gnss_flowgraph.h"
#include "gps_acq_assist.h"


DEFINE_int32(fs_in, 4000000, "Sampling rate, in Ms/s");
DEFINE_int32(max_measurement_duration, 90, "Maximum time waiting for a position fix, in seconds");
DEFINE_int32(num_measurements, 2, "Number of measurements");
DEFINE_string(device_address, "192.168.40.2", "USRP device IP address");
DEFINE_string(subdevice, "A:0", "USRP subdevice");


// For GPS NAVIGATION (L1)
concurrent_queue<Gps_Acq_Assist> global_gps_acq_assist_queue;
concurrent_map<Gps_Acq_Assist> global_gps_acq_assist_map;

std::vector<double> TTFF_v;
const int decimation_factor = 1;

typedef struct  {
    long mtype; // required by SysV message
    double ttff;
} ttff_msgbuf;


class TTFF_GPS_L1_CA_Test: public ::testing::Test
{
public:
    void config_1();
    void config_2();
    void print_TTFF_report(const std::vector<double> & ttff_v, std::shared_ptr<ConfigurationInterface> config_);

    std::shared_ptr<InMemoryConfiguration> config;
    std::shared_ptr<FileConfiguration> config2;

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
};


void TTFF_GPS_L1_CA_Test::config_1()
{
    config = std::make_shared<InMemoryConfiguration>();

    config->set_property("GNSS-SDR.internal_fs_hz", std::to_string(FLAGS_fs_in));

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
    config->set_property("SignalSource.subdevice", FLAGS_subdevice);
    config->set_property("SignalSource.samples", std::to_string(FLAGS_fs_in * FLAGS_max_measurement_duration));
    config->set_property("SignalSource.device_address", FLAGS_device_address);

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
    config->set_property("TelemetryDecoder_1C.decimation_factor", std::to_string(decimation_factor));

    // Set Observables
    config->set_property("Observables.implementation", "GPS_L1_CA_Observables");
    config->set_property("Observables.dump", "false");
    config->set_property("Observables.dump_filename", "./observables.dat");

    // Set PVT
    config->set_property("PVT.implementation", "GPS_L1_CA_PVT");
    config->set_property("PVT.averaging_depth", std::to_string(averaging_depth));
    config->set_property("PVT.flag_averaging", "true");
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


void TTFF_GPS_L1_CA_Test::config_2()
{
    std::string path = std::string(TEST_PATH);
    std::string filename = path + "../../conf/gnss-sdr_GPS_L1_USRP_X300_realtime.conf";
    config2 = std::make_shared<FileConfiguration>(filename);
    config2->set_property("SignalSource.samples", std::to_string(FLAGS_fs_in * FLAGS_max_measurement_duration));
}


void receive_msg()
{
    ttff_msgbuf msg;
    ttff_msgbuf msg_stop;
    msg_stop.mtype = 1;
    msg_stop.ttff = -200.0;
    double ttff_msg = 0.0;
    int msgrcv_size = sizeof(msg.ttff);
    int msqid;
    int msqid_stop = -1;
    key_t key = 1101;
    key_t key_stop = 1102;
    // wait for the queue to be created
    if((msqid = msgget(key, 0644)) == -1){}

    if (msgrcv(msqid, &msg, msgrcv_size, 1, 0) != -1)
        {
            ttff_msg = msg.ttff;
            if( (ttff_msg != 0) && (ttff_msg != -1))
                {
                    TTFF_v.push_back(ttff_msg / (1000.0 / decimation_factor) );
                    LOG(INFO) << "Valid Time-To-First-Fix: " << ttff_msg / (1000.0 / decimation_factor ) << "[s]";
                    // Stop the receiver
                    //while(((msqid_stop = msgget(key_stop, 0644 | IPC_CREAT))) == -1){}
                    while(((msqid_stop = msgget(key_stop, 0644))) == -1){}
                    double msgsend_size = sizeof(msg_stop.ttff);
                    msgsnd(msqid_stop, &msg_stop, msgsend_size, IPC_NOWAIT);
                }

            if(ttff_msg != -1)
                {
                    receive_msg();
                }
        }
    return;
}


void TTFF_GPS_L1_CA_Test::print_TTFF_report(const std::vector<double> & ttff_v, std::shared_ptr<ConfigurationInterface> config_)
{
    std::ofstream ttff_report_file;
    std::string filename = "ttff_report";
    std::string filename_;
    bool time_tag_name = true;

    time_t rawtime;
    struct tm * timeinfo;
    time ( &rawtime );
    timeinfo = localtime ( &rawtime );

    if (time_tag_name)
        {
            std::stringstream strm0;
            const int year = timeinfo->tm_year - 100;
            strm0 << year;
            const int month = timeinfo->tm_mon + 1;
            if(month < 10)
                {
                    strm0 << "0";
                }
            strm0 << month;
            const int day = timeinfo->tm_mday;
            if(day < 10)
                {
                    strm0 << "0";
                }
            strm0 << day << "_";
            const int hour = timeinfo->tm_hour;
            if(hour < 10)
                {
                    strm0 << "0";
                }
            strm0 << hour;
            const int min = timeinfo->tm_min;
            if(min < 10)
                {
                    strm0 << "0";
                }
            strm0 << min;
            const int sec = timeinfo->tm_sec;
            if(sec < 10)
                {
                    strm0 << "0";
                }
            strm0 << sec;

            filename_ = filename + "_" +  strm0.str() + ".txt";
        }
    else
        {
            filename_ = filename + ".txt";
        }

    ttff_report_file.open(filename_.c_str());

    std::vector<double> ttff = ttff_v;
    bool read_ephemeris;
    read_ephemeris = config_->property("GNSS-SDR.SUPL_read_gps_assistance_xml", "false");
    bool agnss;
    agnss = config_->property("GNSS-SDR.SUPL_gps_enabled", "false");
    double sum = std::accumulate(ttff.begin(), ttff.end(), 0.0);
    double mean = sum / ttff.size();
    double sq_sum = std::inner_product(ttff.begin(), ttff.end(), ttff.begin(), 0.0);
    double stdev = std::sqrt(sq_sum / ttff.size() - mean * mean);
    auto max_ttff = std::max_element(std::begin(ttff), std::end(ttff));
    auto min_ttff = std::min_element(std::begin(ttff), std::end(ttff));
    std::string source;
    std::string default_str = "default";
    source = config_->property("SignalSource.implementation", default_str);

    if (ttff_report_file.is_open())
        {
            ttff_report_file << "---------------------------" << std::endl;
            ttff_report_file << " Time-To-First-Fix Report" << std::endl;
            ttff_report_file <<  "---------------------------" << std::endl;
            ttff_report_file << "Initial receiver status: ";
            if (read_ephemeris)
                {
                    ttff_report_file << "Hot start." << std::endl;
                }
            else
                {
                    ttff_report_file << "Cold start." << std::endl;
                }
            ttff_report_file << "A-GNSS: ";
            if (agnss && read_ephemeris)
                {
                    ttff_report_file << "Enabled." << std::endl;
                }
            else
                {
                    ttff_report_file << "Disabled." << std::endl;
                }
            ttff_report_file << "Valid measurements (" << ttff.size() << "/" << FLAGS_num_measurements << "): ";
            for(double ttff_ : ttff) ttff_report_file << ttff_ << " ";
            ttff_report_file << std::endl;
            ttff_report_file << "TTFF mean: " << mean << " [s]" << std::endl;
            if (ttff.size() > 0)
                {
                    ttff_report_file << "TTFF max: " << *max_ttff << " [s]" << std::endl;
                    ttff_report_file << "TTFF min: " << *min_ttff << " [s]" << std::endl;
                }
            ttff_report_file << "TTFF stdev: " << stdev << " [s]" << std::endl;
            ttff_report_file << "Operating System: " << std::string(HOST_SYSTEM) << std::endl;
            ttff_report_file << "Navigation mode: " << "3D" << std::endl;

            if(source.compare("UHD_Signal_Source"))
                {
                    ttff_report_file << "Source: File" << std::endl;
                }
            else
                {
                    ttff_report_file << "Source: Live" << std::endl;
                }
            ttff_report_file << "---------------------------" << std::endl;
        }
    ttff_report_file.close();
    std::cout << "---------------------------" << std::endl;
    std::cout << " Time-To-First-Fix Report" << std::endl;
    std::cout << "---------------------------" << std::endl;
    std::cout << "Initial receiver status: ";
    if (read_ephemeris)
        {
            std::cout << "Hot start." << std::endl;
        }
    else
        {
            std::cout << "Cold start." << std::endl;
        }
    std::cout << "A-GNSS: ";
    if (agnss && read_ephemeris)
        {
            std::cout << "Enabled." << std::endl;
        }
    else
        {
            std::cout << "Disabled." << std::endl;
        }
    std::cout << "Valid measurements (" << ttff.size() << "/" << FLAGS_num_measurements << "): ";
    for(double ttff_ : ttff) std::cout << ttff_ << " ";
    std::cout << std::endl;
    std::cout << "TTFF mean: " << mean << " [s]" << std::endl;
    if (ttff.size() > 0)
        {
            std::cout << "TTFF max: " << *max_ttff << " [s]" << std::endl;
            std::cout << "TTFF min: " << *min_ttff << " [s]" << std::endl;
        }
    std::cout << "TTFF stdev: " << stdev << " [s]" << std::endl;
    std::cout << "Operating System: " << std::string(HOST_SYSTEM) << std::endl;
    std::cout << "Navigation mode: " << "3D" << std::endl;

    if(source.compare("UHD_Signal_Source"))
        {
            std::cout << "Source: File" << std::endl;
        }
    else
        {
            std::cout << "Source: Live" << std::endl;
        }
    std::cout << "---------------------------" << std::endl;
}


TEST_F(TTFF_GPS_L1_CA_Test, ColdStart)
{
    unsigned int num_measurements = 0;

    config_1();
    // Ensure Cold Start
    config->set_property("GNSS-SDR.SUPL_gps_enabled", "false");
    config->set_property("GNSS-SDR.SUPL_read_gps_assistance_xml", "false");

    config_2();
    // Ensure Cold Start
    config2->set_property("GNSS-SDR.SUPL_gps_enabled", "false");
    config2->set_property("GNSS-SDR.SUPL_read_gps_assistance_xml", "false");
    config2->set_property("PVT.flag_rtcm_server", "false");

    for(int n = 0; n < FLAGS_num_measurements; n++)
        {
            // Create a new ControlThread object with a smart pointer
            std::unique_ptr<ControlThread> control_thread(new ControlThread(config));

            // record startup time
            struct timeval tv;
            gettimeofday(&tv, NULL);
            long long int begin = tv.tv_sec * 1000000 + tv.tv_usec;

            std::cout << "Starting measurement " << num_measurements + 1 << " / " << FLAGS_num_measurements << std::endl;

            // start receiver
            try
            {
                    control_thread->run();
            }
            catch( boost::exception & e )
            {
                    std::cout << "Boost exception: " << boost::diagnostic_information(e);
            }
            catch(std::exception const&  ex)
            {
                    std::cout  << "STD exception: " << ex.what();
            }

            // stop clock
            gettimeofday(&tv, NULL);
            long long int end = tv.tv_sec * 1000000 + tv.tv_usec;
            double ttff = static_cast<double>(end - begin) / 1000000.0;

            std::shared_ptr<GNSSFlowgraph> flowgraph = control_thread->flowgraph();
            EXPECT_FALSE(flowgraph->running());

            num_measurements = num_measurements + 1;
            std::cout << "Just finished measurement " << num_measurements << ", which took " << ttff << " seconds." << std::endl;
            if(n < FLAGS_num_measurements - 1)
                {
                    std::srand(std::time(0)); // use current time as seed for random generator
                    int random_variable = std::rand();
                    float random_variable_0_1 = static_cast<float>(random_variable) / static_cast<float>( RAND_MAX );
                    int random_delay_s = static_cast<int>(random_variable_0_1 * 25.0);
                    std::cout << "Waiting a random amount of time (from 5 to 30 s) to start a new measurement... " << std::endl;
                    std::cout << "This time will wait " << random_delay_s + 5 << " s." << std::endl << std::endl;
                    std::this_thread::sleep_until(std::chrono::system_clock::now() + std::chrono::seconds(5) + std::chrono::seconds(random_delay_s));
                }
        }

    // Print TTFF report
    print_TTFF_report(TTFF_v, config);
}


TEST_F(TTFF_GPS_L1_CA_Test, HotStart)
{
    unsigned int num_measurements = 0;
    TTFF_v.clear();

    config_1();
    // Ensure Hot Start
    config->set_property("GNSS-SDR.SUPL_gps_enabled", "true");
    config->set_property("GNSS-SDR.SUPL_read_gps_assistance_xml", "true");

    for(int n = 0; n < FLAGS_num_measurements; n++)
        {
            // Create a new ControlThread object with a smart pointer
            std::unique_ptr<ControlThread> control_thread(new ControlThread(config));

            // record startup time
            struct timeval tv;
            gettimeofday(&tv, NULL);
            long long int begin = tv.tv_sec * 1000000 + tv.tv_usec;

            std::cout << "Starting measurement " << num_measurements + 1 << " / " << FLAGS_num_measurements << std::endl;

            // start receiver
            try
            {
                    control_thread->run();
            }
            catch( boost::exception & e )
            {
                    std::cout << "Boost exception: " << boost::diagnostic_information(e);
            }
            catch(std::exception const&  ex)
            {
                    std::cout  << "STD exception: " << ex.what();
            }

            // stop clock
            gettimeofday(&tv, NULL);
            long long int end = tv.tv_sec * 1000000 + tv.tv_usec;
            double ttff = static_cast<double>(end - begin) / 1000000.0;

            std::shared_ptr<GNSSFlowgraph> flowgraph = control_thread->flowgraph();
            EXPECT_FALSE(flowgraph->running());

            num_measurements = num_measurements + 1;
            std::cout << "Just finished measurement " << num_measurements << ", which took " << ttff << " seconds." << std::endl;
            if(n < FLAGS_num_measurements - 1)
                {
                    std::srand(std::time(0)); // use current time as seed for random generator
                    int random_variable = std::rand();
                    float random_variable_0_1 = static_cast<float>(random_variable) / static_cast<float>( RAND_MAX );
                    int random_delay_s = static_cast<int>(random_variable_0_1 * 25.0);
                    std::cout << "Waiting a random amount of time (from 5 to 30 s) to start new measurement... " << std::endl;
                    std::cout << "This time will wait " << random_delay_s + 5 << " s." << std::endl << std::endl;
                    std::this_thread::sleep_until(std::chrono::system_clock::now() + std::chrono::seconds(5) + std::chrono::seconds(random_delay_s));
                }
        }

    // Print TTFF report
    print_TTFF_report(TTFF_v, config);
}


int main(int argc, char **argv)
{
    std::cout << "Running Time-To-First-Fix test..." << std::endl;
    int res = 0;
    TTFF_v.clear();
    testing::InitGoogleTest(&argc, argv);
    google::ParseCommandLineFlags(&argc, &argv, true);
    google::InitGoogleLogging(argv[0]);

    // Create SysV message queue to read TFFF measurements
    key_t sysv_msg_key;
    int sysv_msqid;
    sysv_msg_key = 1101;
    int msgflg = IPC_CREAT | 0666;
    if ((sysv_msqid = msgget(sysv_msg_key, msgflg )) == -1)
    {
        std::cout<<"GNSS-SDR can not create message queues!" << std::endl;
        throw new std::exception();
    }

    // Start queue thread
    std::thread receive_msg_thread(receive_msg);

    // Run the Tests
    try
    {
            res = RUN_ALL_TESTS();
    }
    catch(...)
    {
            LOG(WARNING) << "Unexpected catch";
    }

    // Terminate the queue thread
    ttff_msgbuf msg;
    msg.mtype = 1;
    msg.ttff = -1;
    int msgsend_size;
    msgsend_size = sizeof(msg.ttff);
    msgsnd(sysv_msqid, &msg, msgsend_size, IPC_NOWAIT);
    receive_msg_thread.join();
    msgctl(sysv_msqid, IPC_RMID, NULL);
    google::ShutDownCommandLineFlags();
    return res;
}
