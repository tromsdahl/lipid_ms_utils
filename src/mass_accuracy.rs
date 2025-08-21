use std::io;

pub fn mass_accuracy() {
    loop {
        println!("\n");
        println!("Enter in the measured mass vs. calculated mass to determine the ppm mass error:");
        println!("Tip: Enter in the measured mass first then calculated mass separated by /");

        let mut mass_entry = String::new();

        io::stdin()
            .read_line(&mut mass_entry)
            .expect("Failed to read line.");

        let measured_index = match mass_entry.find('/') {
            Some(n) => n,
            None => {
                println!("Formatting not recognized! A '/' wasn't found!");
                continue;
            }
        };
        let measured_mass: f64 = match mass_entry[..measured_index].parse() {
            Ok(n) => n,
            Err(e) => {
                println!("The measured mass wasn't recognized! You entered: {e}");
                continue;
            }
        };

        if !(measured_mass > 0.0) {
            println!("The measured mass is not a positive number! You entered: {measured_mass}");
            continue;
        }

        let calculated_index = mass_entry.len();
        let calculated_mass: f64 = match mass_entry[measured_index + 1..calculated_index].trim().parse() {
            Ok(n) => n,
            Err(e) => {
                println!("The measured mass wasn't recognized! You entered: {e}");
                continue;
            }
        };

        if !(calculated_mass > 0.0) {
            println!("The expected mass is not a positive number! You entered: {calculated_mass}");
            continue;
        }

        println!("\n");
        println!("Measured mass: {measured_mass}");
        println!("Calculated mass: {calculated_mass}");

        let ppm_error = ((measured_mass - calculated_mass) / calculated_mass) * 1e6;

        println!("\n");
        println!("Mass error: {ppm_error:.2} ppm");

        break;
    };

}