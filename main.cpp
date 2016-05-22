//
//  main.cpp
//  pata-risk
//
//  Created by a1091793 on 5/05/2016.
//  Copyright © 2016 University of Adelaide and Bushfire and Natural Hazards CRC. All rights reserved.
//

#include <iostream>
#include <sstream>
#include <string>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <blink/raster/utility.h>
#include <blink/iterator/zip_range.h>
#include "CCCALM.hpp"
#include "InundateLandscape.h"
#include "ReadInControlsAndGuages.h"
#include "ReadGraphsFromFile.h"
#include "CalcDamage.hpp"
#include "AggregateMap.h"

int main(int argc, const char * argv[]) {
    // insert code here...
    namespace prog_opt = boost::program_options;
    namespace raster_it = blink::iterator;
    namespace raster_util = blink::raster;
    
    
    //	    namespace raster_util = blink::raster;
    
    namespace fs = boost::filesystem;
    fs::path run_path = fs::current_path();
    
    /**********************************/
    /*        Program options         */
    /**********************************/
    
    
    ////VARIABLES FOR LANDUSE MODELLING
    int device_id = -1;
    std::string project_path = "current dir";
    std::string parameter_path = "parameters.txt";
    int numSteps = 10;
    
    
    
    /////VARIABLES FOR INUNDATION MODELLING
    // Need to specify elevation grid
    // Need to specify channel
    //	std::string feature_file;
    std::string dem_file;
    std::string sclimate_sc("no_Sc");
    //	std::string changes_file;
    //	std::string log_file;
    //    std::string fd_file;
    std::string hydro_paths_file;
    std::string output_file;
    std::string channel_graph_file;
    std::string controls_file("no_file");
    std::string guage_stem;
    std::string mask_file;
    //	bool do_print = false;
    //	std::string file_print;
    //    unsigned int trim_level;
    bool do_interp_hgl = false;
    bool write_controls= false;
    double outflow_levels = 1.0;
    double source_levels = 0.0;
    bool use_controls_file = false;
    double rate = 0.06;
    std::string outfile;
    
    
    
    prog_opt::options_description desc("Allowed options");
    desc.add_options()
    ("help,h", "produce help message")
    //////////////    LANDUSE     //////////////
    ("project-path,j", prog_opt::value<std::string>(&(project_path))->default_value("no_dir"), "path of the directory with required files")
    //        ("flow-dir-map,i", prog_opt::value<std::string>(&fd_file), "path of the gdal capatible flow direction data file")
    ("parameters,p", prog_opt::value<std::string>(&(parameter_path))->default_value("no_file"), "the path of parameters file")
    ("timesteps,n", prog_opt::value<int>(&(numSteps))->default_value(10), "the number of CA steps to perform.")
    
    //////////////    INUNDATION     //////////////
    ("dem-map,d", prog_opt::value<std::string>(&dem_file), "path of the gdal capatible elevation data file")
    //        ("flow-dir-map,i", prog_opt::value<std::string>(&fd_file), "path of the gdal capatible flow direction data file")
    ("hydro-paths-file,t", prog_opt::value<std::string>(&hydro_paths_file)->default_value("hydro-paths.tif"), "path of the output map where each pixel is assigned the location on channel that the pixel is hydrologically connected to ")
    ("mask,m", prog_opt::value<std::string>(&mask_file)->default_value("mask.tiff"), "Path to the mask file - no losses where value = 0 in this mask as this is where the water should be -= e.g. in the channel etc")
    ("channel-graph,g", prog_opt::value<std::string>(&channel_graph_file), "path of the graphml representation of the channel")
    ("climate-sc,c", prog_opt::value<std::string>(&sclimate_sc), "path of text file with guage location and levels, takes precedence over channel graph")
    ("guages-stem,a", prog_opt::value<std::string>(&guage_stem), "path of the directory where the guages files are located")
//    ("output-flood-height-file,o", prog_opt::value<std::string>(&sclimate_sc)->default_value("rcp45"), "path of the output map where each pixel is assigned the flood height at that pixel")
    //        ("trim-branches,t", prog_opt::value<unsigned int>(&trim_level)->default_value(0), "remove branches if they are composed with a number of pixels less than this amount")
//    ("read-controls,r", prog_opt::value<std::string>(&controls_file), "Name of controls file to read into. Take precedence over guage table file")
    ("outflow-levels,u", prog_opt::value<double>(&outflow_levels), "Flood level assigned to outflow nodes (least precedence)")
    ("source-levels,s", prog_opt::value<double>(&source_levels), "Flood level assigned to outflow nodes (least precedence)")
    ("rate,r", prog_opt::value<double>(&rate)->default_value(0.06), "Discount rate used for present value analysis")
    ("out,o", prog_opt::value<std::string>(&outfile)->default_value("pv_risk.txt"), "file for printing the present value of risk to");
    
    
    
    prog_opt::variables_map vm;
    prog_opt::store(prog_opt::parse_command_line(argc, argv, desc), vm);
    prog_opt::notify(vm);
    if (vm.count("help"))
    {
        std::cout << desc << "\n";
        return 1;
    }
//    
//    std::cout << "************************************" << std::endl;
//    std::cout << "Constrained CA for land-use dynamics" << std::endl;
    
    initialiseCCCALM(project_path, parameter_path, numSteps, device_id );
    setVerbose();
    
    
    
//    std::cout << "End of land-use simulation" << std::endl;
    //    CCCALM::closeCA(data);
    //    system( "pause" );
    
    fs::path channel_graph_path(channel_graph_file);
    fs::path dem_file_path(dem_file);
    //	fs::path changes_file_path(changes_file);
    fs::path guage_table_path(guage_stem);
    //  fs::path fd_file_path(fd_file);
    fs::path hydro_paths_file_path(hydro_paths_file);
    fs::path output_file_path(output_file);
    fs::path controls_file_path(controls_file);
    fs::path guages_stem_path(guage_stem);
    fs::path mask_path(mask_file);
    
    // Check file exists
    if (!fs::exists(channel_graph_path))
    {
        std::stringstream ss;
        ss << channel_graph_path << " does not exist";
        throw std::runtime_error(ss.str());
        return (EXIT_FAILURE);
    }
    
    if (!fs::exists(dem_file_path))
    {
        std::stringstream ss;
        ss << dem_file_path << " does not exist";
        throw std::runtime_error(ss.str());
        return (EXIT_FAILURE);
    }
    
    if (!fs::exists(hydro_paths_file_path))
    {
        std::stringstream ss;
        ss << hydro_paths_file_path << " does not exist";
        throw std::runtime_error(ss.str());
        return (EXIT_FAILURE);
    }
    
    if (!fs::exists(guages_stem_path))
    {
        std::stringstream ss;
        ss << guages_stem_path << " does not exist";
        throw std::runtime_error(ss.str());
        return (EXIT_FAILURE);
    }
    
    if (!fs::exists(mask_path))
    {
        std::stringstream ss;
        ss << mask_path << " does not exist";
        throw std::runtime_error(ss.str());
        return (EXIT_FAILURE);
    }
    
    
    ControlsSPtr controls;
    if(use_controls_file)
    {
        if (!fs::exists(controls_file_path))
        {
            std::stringstream ss;
            ss << controls_file_path << " does not exist";
            throw std::runtime_error(ss.str());
            return (EXIT_FAILURE);
        }
        controls = readInControls(controls_file_path);
    }
    
    
    /**********************************/
    /*       Create graph object      */
    /**********************************/
    Graph channel_grph;
    
    
    /**********************************/
    /*         Read in Graph           */
    /**********************************/
    std::cout << "\n\n*************************************\n";
    std::cout <<     "*             Read in Graphs          *\n";
    std::cout <<     "*************************************" << std::endl;
    //    readGraphFromFile(control_graph_path, control_grph);
    readGraphFromFile(channel_graph_path, channel_grph);
    
    auto dem = raster_util::open_gdal_raster<double>(dem_file_path.string(), GA_ReadOnly);
    auto hydro_connect = raster_util::open_gdal_raster<int>(hydro_paths_file_path.string(), GA_ReadOnly);
//    auto inundation = raster_util::create_gdal_raster_from_model<double>(output_file, dem);
    
    auto mask_map = raster_util::open_gdal_raster<int>(mask_path, GA_ReadOnly);
    
    
    
    std::vector<int> aris = {20, 50, 100, 200, 500, 1000};
    std::vector<double> freqs(aris.size());
    for (int i = 0; i < aris.size(); ++i) {
        freqs[i] = 1 / double(aris[i]);
    }
    
    std::vector<double> risk_by_year;
    
    fs::path pv_risk_map_file = run_path / ("pv_risk.tif");
    DoubleRaster pv_risk_raster;

    for (int step = 0; step < numSteps; ++step)
    {
        IntRaster& landuse_map = stepCCCALM(1);
        std::vector<DoubleRaster> loss_maps(6);
        
        if (step == 0)
        {
          pv_risk_raster  = raster_util::create_gdal_raster_from_model<double>(pv_risk_map_file, dem);
        }
//        
//        
//        // Calculate Hazard, generating hazard map
        for(int ari = 0; ari < aris.size(); ++ari)
        {
//
            GuagesSPtr guages;
            guage_table_path = guages_stem_path / (std::to_string(2016+step) + "_" + sclimate_sc + "_ari" + std::to_string(aris[ari]) + ".txt");
            if (!fs::exists(guage_table_path))
            {
                std::stringstream ss;
                ss << guage_table_path << " does not exist";
                throw std::runtime_error(ss.str());
                return (EXIT_FAILURE);
            }
            guages = readInGuages(guage_table_path);
//
            fs::path output_file_path = run_path / ("inundation_" + std::to_string(2016+step) + "_" + sclimate_sc + "_ari" + std::to_string(aris[ari]) + ".tif");
            fs::path output_file_path_depth = run_path / ("inundation_depth_" + std::to_string(2016+step) + "_" + sclimate_sc + "_ari" + std::to_string(aris[ari]) + ".tif");
            fs::path output_file_path_prop = run_path / ("inundation_prop_" + std::to_string(2016+step) + "_" + sclimate_sc + "_ari" + std::to_string(aris[ari]) + ".tif");
            fs::path loss_map_path = run_path / ("loss_" + std::to_string(2016+step) + "_" + sclimate_sc + "_ari" + std::to_string(aris[ari]) + ".tif");
            auto inundation = raster_util::create_gdal_raster_from_model<double>(output_file_path, dem);
            const_cast<GDALRasterBand *>(inundation.get_gdal_band())->SetNoDataValue(0.0);
            inundateLandscape(inundation, dem, hydro_connect, channel_grph, guages, controls);
            auto agg_map_depth = raster_util::create_gdal_raster_from_model<double>(output_file_path_depth, landuse_map);
            auto agg_map_proportion = raster_util::create_gdal_raster_from_model<double>(output_file_path_prop, landuse_map);
            aggregateMaps(inundation, landuse_map, agg_map_depth, agg_map_proportion);
            loss_maps[ari] = raster_util::create_gdal_raster_from_model<double>(loss_map_path, landuse_map);
            calcNetLosses(agg_map_depth, agg_map_proportion, mask_map, landuse_map, loss_maps[ari]);
        }
    
        
        // create a raster data set, with same dimensions as map1
        fs::path risk_raster_file = run_path / ("risk_" + std::to_string(step) + ".tif");
        DoubleRaster risk_raster = raster_util::create_gdal_raster_from_model<double>(risk_raster_file, dem);
        
        
        
//        risk_raster.setNoDataValue(0.0);
        double net_risk = 0;
        
        double risk_i;
        
        auto zip = raster_it::make_zip_range(std::ref(loss_maps[0]), std::ref(loss_maps[1]), std::ref(loss_maps[2]), std::ref(loss_maps[3]), std::ref(loss_maps[4]), std::ref(loss_maps[5]), std::ref(risk_raster), std::ref(pv_risk_raster));
        for (auto i : zip)
        {
            const double & loss1 = std::get<0>(i);
            const double & loss2 = std::get<1>(i);
            const double & loss3 = std::get<2>(i);
            const double & loss4 = std::get<3>(i);
            const double & loss5 = std::get<4>(i);
            const double & loss6 = std::get<5>(i);
            auto & risk = std::get<6>(i);
            auto & pv_risk = std::get<7>(i);
            
            const double exp_loss1 = freqs[0] * loss1;
            const double exp_loss2 = freqs[1] * loss2;
            const double exp_loss3 = freqs[2] * loss3;
            const double exp_loss4 = freqs[3] * loss4;
            const double exp_loss5 = freqs[4] * loss5;
            const double exp_loss6 = freqs[5] * loss6;
            
//            risk_i = ( (0     + exp_loss1) * std::abs(1       - freq[0])
//                       (exp_loss1 + exp_loss2) * std::abs(freq[0] - freq[1]) +
//                       (exp_loss2 + exp_loss3) * std::abs(freq[1] - freq[2]) +
//                       (exp_loss3 + exp_loss4) * std::abs(freq[2] - freq[3]) +
//                       (exp_loss4 + exp_loss5) * std::abs(freq[3] - freq[4]) +
//                       (exp_loss5 + exp_loss6) * std::abs(freq[4] - freq[5]) +
//                     ) / 2;
            
            
            risk_i = ( (0     + exp_loss1)    * std::abs(loss1 - 0    ) +
                      (exp_loss1 + exp_loss2) * std::abs(loss2 - loss1) +
                      (exp_loss2 + exp_loss3) * std::abs(loss3 - loss2) +
                      (exp_loss3 + exp_loss4) * std::abs(loss4 - loss3) +
                      (exp_loss4 + exp_loss5) * std::abs(loss5 - loss4) +
                      (exp_loss5 + exp_loss6) * std::abs(loss6 - loss5)
                      ) / 2;
            
            
            risk = risk_i;
            net_risk += risk_i;
            pv_risk += risk_i;
        }
        
        risk_by_year.push_back(net_risk);
        
    }
    
    
    // Calc. NPV risk.
    double pv_risk = 0;
    for (int i = 0; i < risk_by_year.size(); ++i)
    {
        pv_risk += risk_by_year[i] / std::pow((1+rate), i);
    }
    
    std::ofstream ofs(outfile.c_str());
    if (ofs.is_open())
    {
        ofs << pv_risk;
        ofs.close();
    }
    
    
    return pv_risk;
}
