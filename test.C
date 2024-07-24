{
    // Including the necessary headers
    #include <vector>
    #include <iostream>

    std::vector<double> bkg = {1.0, 2.0, 3.0}; // Example initialization
    double liveTimeTh = 10.0;
    double liveTimeBkg = 5.0;

    // Check if liveTimeBkg is not zero to avoid division by zero
    if (liveTimeBkg != 0) {
        double scale = liveTimeTh / liveTimeBkg;
        for (auto& value : bkg) {
            value *= scale;
        }
    } else {
        std::cerr << "Error: liveTimeBkg is zero, division by zero is not allowed." << std::endl;
    }

    // Output the rescaled bkg vector for verification
    for (const auto& value : bkg) {
        std::cout << value << " ";
    }
    std::cout << std::endl;
}
