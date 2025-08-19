use std::io;

fn main() {

    println!("Lipid MS Utilities");
    println!("Enter a number to select an option:");

    loop {        
        
        println!("\n");
        println!("1. Lipid mass calculator");
        println!("2. quit.");

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
        1 => crate::lipid_calculator::lipid_calc_funct(),
        2 => {
            println!("Program closing!");
            break;
        },
        _ => continue,
    }

    }

}

mod lipid_calculator {
    use std::io;

    // set atom masses to constants
    const CARBON: f64 = 12.000000;
    const HYDROGEN: f64 = 1.007825;
    const OXYGEN: f64 = 15.994915;
    //const PHOSPH: f64 = 30.973762;
    //const NITROGEN: f64 = 14.003074;
    const AMMONIUM: f64 = 18.03382555;
    const SODIUM: f64 = 22.989769;
    const POTASSIUM: f64 = 38.963707;
    //const SULFUR: f64 = 31.972072;

    pub fn lipid_calc_funct () {
        loop {
            
            let lipid_classes = vec!["TG", "CE", "MG", "DG", "PC", "LPC", "PE", "LPE", "PS", "LPS", "PI", "LPI"];
            
            println!("\n");
            println!("Enter in the lipid species to calculate its accurate mass and adducts (e.g. PC-34:2):");
            println!("Options include: TG, CE, MG, DG, PC, LPC, PE, LPE, PS, LPS, PI, LPI");
            println!("Tip: include a hyphen between the class and total fatty acid carbons and a colon between fatty acid carbons and unsaturations");

            let mut lipid_entry = String::new();

            io::stdin()
                .read_line(&mut lipid_entry)
                .expect("Failed to read line.");

            let class_index = lipid_entry.find('-').unwrap();
            let lipid_class = &lipid_entry[..class_index];

            if !lipid_classes.contains(&lipid_class) {
                println!("The lipid class entered is not recognized! You entered: {lipid_class}");
                continue;
            }

            let carbon_index = lipid_entry.find(':').unwrap();
            let fa_carbon = &lipid_entry[class_index + 1usize..carbon_index];
            let fa_carbon_no: f64 = fa_carbon.parse().unwrap();

            if !(fa_carbon_no > 0.0) {
                println!("The FA carbons is not a number! You entered: {fa_carbon}");
                continue;
            }

            let desat_index = lipid_entry.len();
            let fa_desat = &lipid_entry[carbon_index + 1..desat_index];
            let fa_desat_no: f64 = fa_desat.trim().parse().unwrap();

            if !(fa_desat_no >= 0.0) {
                println!("The FA unsaturations is not a number! You entered: {fa_desat}");
                continue;
            }

            println!("\n");
            println!("Lipid class: {lipid_class}");
            println!("Total no. of carbons: {fa_carbon}");
            println!("Total no. of unsaturations: {fa_desat}");

            let calc_masses = match lipid_class {
                "TG" => tg_calc(&fa_carbon_no, &fa_desat_no),
                _ => {
                    println!("That lipid class hasn't been defined yet! Hold tight buckaroo!");
                    continue;
                },
            
            };

            let monoiso_mass = &calc_masses[0];
            let h_mass = &calc_masses[1];
            let na_mass = &calc_masses[2];
            let k_mass = &calc_masses[3];
            let nh4_mass = &calc_masses[4];

            println!("\n");
            println!("Monoisotopic mass: {monoiso_mass}");
            println!("[M+H]+: {h_mass}");
            println!("[M+Na]+: {na_mass}");
            println!("[M+K]+: {k_mass}");
            println!("[M+NH4]+: {nh4_mass}");

            break;
        };

            fn tg_calc(carb: &f64, desat: &f64) -> Vec<String> {
                let mut calculated_masses = Vec::new();

                let monoiso_mass: f64 = (CARBON * (carb + 3.0)) + (HYDROGEN * ((desat * 2.0) - (desat * 2.0) + 2.0)) + (OXYGEN * 6.0);

                let tg_h_adduct = &monoiso_mass + HYDROGEN;
                let tg_na_adduct = &monoiso_mass + SODIUM;
                let tg_k_adduct = &monoiso_mass + POTASSIUM;
                let tg_nh4_adduct = &monoiso_mass + AMMONIUM;

                let monoiso_mass_form = format!("{monoiso_mass:.4}");
                let h_adduct_form = format!("{tg_h_adduct:.4}");
                let na_adduct_form = format!("{tg_na_adduct:.4}");
                let k_adduct_form = format!("{tg_k_adduct:.4}");
                let nh4_adduct_form = format!("{tg_nh4_adduct:.4}");

                calculated_masses.push(monoiso_mass_form);
                calculated_masses.push(h_adduct_form);
                calculated_masses.push(na_adduct_form);
                calculated_masses.push(k_adduct_form);
                calculated_masses.push(nh4_adduct_form);

                calculated_masses
            }
    }
}
