use std::io;

pub fn molec_mass_calc() {
    loop {
        println!("\n");
        println!("Enter in a molecular formula to calculate its accurate monoisotopic mass.");

        let mut formula = String::new();

        io::stdin()
            .read_line(&mut formula)
            .expect("Failed to read line.");

        let calculated_mass = match calculate_mass(&formula.trim()) {
            Ok(m) => m,
            Err(e) => {
                println!("There was an error! {e}");
                continue;
            }
        };

        let protonated_adduct = 1.007825 + &calculated_mass;
        let sodiated_adduct = 22.989769 + &calculated_mass;
        let minus_proton = &calculated_mass - 1.007825;

        println!("\n");
        println!("Monoisotopic mass is: {calculated_mass:.4}");
        println!("\n");
        println!("[M+H]: {protonated_adduct:.4}");
        println!("[M+Na]: {sodiated_adduct:.4}");
        println!("[M-H]: {minus_proton:.4}");
        break;
}

        fn calculate_mass(formula: &str) -> Result<f64, String> {
            // use a counter to keep track of the mass as you're reading through the molecular formula
            let mut total_mass = 0.0;

            // use peekable iterator to see next char without consuming
            let mut chars = formula.chars().peekable();

            // loop through formula to consume all chars
            while let Some(ch) = chars.next() {
                
                // first check that the first letter is a capital letter
                if !ch.is_ascii_uppercase() {
                    return Err(
                        format!("Formatting not correct! There should be an uppercase letter at the beginning. You entered:{ch}")
                    );
                }

                // check next char to see if number, atom, or second (lowercase) letter to atom
                let mut element_symbol = ch.to_string();
                if let Some(&next_ch) = chars.peek() { // note that this only executes if the next letter is lowercase, otherwise the element_symbol remains the same the code after this section runs
                    if next_ch.is_ascii_lowercase() {
                        // add lowercase to original letter for complete atom symbol
                        element_symbol.push(chars.next().unwrap());
                    }
                }

                // find element in list of available elements from another function (element_lookup)
                let mass = match element_lookup(&element_symbol) {
                    Some(m) => m,
                    None => return Err(
                        format!("Unknown element symbol:{element_symbol}")
                    ),
                };

                // now see how many of that element there are
                let mut count_str = String::new();

                // will add each number to string (e.g. 3 then 2 for 32)
                while let Some(&next_ch) = chars.peek() {
                    if next_ch.is_ascii_digit() {
                        count_str.push(chars.next().unwrap());
                    } else {
                        break;
                    }
                }

                let count: u32 = if count_str.is_empty() {
                    1
                } else {
                    match count_str.parse() {
                        Ok(n) => n,
                        Err(_) => return Err(
                            format!("Invalid number format for count: {count_str}")
                        ),
                    }
                };

                // add the element's mass to the runing total from above
                total_mass += mass * (count as f64);


            };

            Ok(total_mass)

        }

        fn element_lookup(symbol: &str) -> Option<f64> {
            match symbol {
                // organic atoms
                "C" => Some(12.000000),
                "H" => Some(1.007825),
                "N" => Some(14.003074),
                "O" => Some(15.994915),
                "P" => Some(30.973762),
                "S" => Some(31.972071),

                // halogens
                "F" => Some(18.998403),
                "Cl" => Some(34.968853),
                "Br" => Some(78.918338),
                "I" => Some(126.904473),

                // metals
                "Na" => Some(22.989769),
                "K" => Some(38.963706),
                "Mg" => Some(23.985042),
                "Ca" => Some(39.962591),
                _=> None // if element is not found or hasn't been added yet
            }
        }
}