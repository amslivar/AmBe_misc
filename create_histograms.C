//to run, do
//root -l -q 'create_histograms.C("AmBe_Cal_Source_Cal_with_Atten_15mins.dat", "AmBe_Run_with_atten_6hrs.dat", "AmBe_Run_without_atten_7hrs.dat", "backgroundrunwithoutatten_04182024.dat", "background_run_with_atten_04232024.dat", "Th228_Cal.dat", "Th228_Cal_with_Atten.dat", "Th228_Cal_without_Atten.dat")'
//("AmBe_Cal_Source_Cal_with_Atten_15mins.dat", "AmBe_Run_with_atten_6hrs.dat", "AmBe_Run_without_atten_7hrs.dat", "backgroundrunwithoutatten_04182024.dat", "background_run_with_atten_04232024.dat", "Th228_Cal.dat", "Th228_Cal_with_Atten.dat", "Th228_Cal_without_Atten.dat")
  // Add your .dat files here
#include <TMath.h>
#include <TF1.h>

void create_histograms() {

    //std::vector<std::string> datFiles = {"AmBe_Cal_Source_Cal_with_Atten_15mins.dat", "AmBe_Run_with_atten_6hrs.dat", "AmBe_Run_without_atten_7hrs.dat", "backgroundrunwithoutatten_04182024.dat", "background_run_with_atten_04232024.dat", "Th228_Cal.dat", "Th228_Cal_with_Atten.dat", "Th228_Cal_without_Atten.dat"};
    std::vector<std::string> datFiles = {"AmBe_Cal_Source_Cal_with_Atten_15mins.dat"}; 


    // Create a canvas to draw histograms
    TCanvas* canvas = new TCanvas("canvas", "Histograms", 800, 600);
    
    // Create a PDF to save all histograms
    canvas->Print("histograms.pdf[");

    for (const auto& fileName : datFiles) {
        // Open the .dat file
        ifstream infile(fileName.c_str());
        if (!infile) {
            std::cerr << "Error opening file: " << fileName << std::endl;
            continue;
        }

        // Create a histogram
        TH1F* hist = new TH1F(fileName.c_str(), fileName.c_str(), 100, 0, 6000);  // Customize bins and range as needed

        double value;
        while (infile >> value) {
            hist->Fill(value);
        }
        infile.close();

        // Draw the histogram on the canvas
        hist->Draw();

        // Add the histogram to the PDF
        canvas->Print("histograms.pdf");
    }

    // Close the PDF file
    canvas->Print("histograms.pdf]");
}

void run_create_histograms_from_dat_files(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage: root -l -q 'macro.C(\"file1.dat\", \"file2.dat\", ...)' " << std::endl;
        return;
    }
    
    std::vector<std::string> datFiles;
    for (int i = 1; i < argc; ++i) {
        datFiles.push_back(argv[i]);
    }

    create_histograms();
}
