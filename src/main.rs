use std::io::{self, Read};

fn main() {

    println!("Lipid MS Utilities");
    println!("Enter a number to select an option:");

    loop {        
        
        println!("1. Lipid mass calculator");
        println!("Press 'q' to quit.");

        let mut user_selection = String::new();

        io::stdin()
            .read_line(&mut user_selection)
            .expect("Failed to read line.");

        let user_selection: u32 = match user_selection.trim().parse() {
            Ok(num) => num,
            Err(_) => continue,
        };

        println!("You entered: {user_selection}");

    match user_selection {
        1 => lipid_calculator,
        'q' | 'Q' => {
            println!("Program closing!");
            break;
        },
        _ => break,
    }

    }

}

mod lipid_calculator {
    // set atom masses to constants
    const CARBON: f64 = 12.000000;
    const HYDROGEN: f64 = 1.007825;
    const OXYGEN: f64 = 15.994915;
    const PHOSPH: f64 = 30.973762;
    const NITROGEN: f64 = 14.003074;
    const AMMONIUM: f64 = 18.03382555;
    const SODIUM: f64 = 22.989769;
    const POTASSIUM: f64 = 38.963707;
    const SULFUR: f64 = 31.972072;

    fn lipid_calc_funct () {
        loop (
            println!("Enter in the lipid class:");
            println!("Options:");
            println!("TG");

            let mut lipid_class = String::new();

            io::stdin()
                .read_line(&mut lipid_class)
                .expect("Failed to read line.");

            


        )
    }

    fn tag_mw(fa_c: f64, fa_des: f64) -> String {
        let monoiso_mass = (CARBON * (fa_c + 3.0)) + (HYDROGEN * ((fa_c * 2.0) - (fa_des * 2.0) + 2.0)) + (OXYGEN * 6.0);
        let tag_h_adduct = &monoiso_mass + HYDROGEN;
        let tag_na_adduct = &monoiso_mass + SODIUM;
        let tag_k_adduct = &monoiso_mass + POTASSIUM;
        let tag_nh4_adduct = &monoiso_mass + AMMONIUM;

        format!("Monoisotopic mass: {monoiso_mass:.4}")

    }
}
