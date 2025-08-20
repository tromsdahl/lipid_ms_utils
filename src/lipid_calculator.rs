use std::io;

// set atom masses to constants
const CARBON: f64 = 12.000000;
const HYDROGEN: f64 = 1.007825;
const OXYGEN: f64 = 15.994915;
const PHOSPH: f64 = 30.973762;
const NITROGEN: f64 = 14.003074;
const AMMONIUM: f64 = 18.03382555;
const SODIUM: f64 = 22.989769;
const POTASSIUM: f64 = 38.963707;
//const SULFUR: f64 = 31.972072;

pub fn lipid_calc_funct () {
    loop {
        
        let lipid_classes = vec!["TG", "CE", "MG", "DG", "PC", "LPC", "PE", "LPE", "PS", "LPS", "PI", "LPI", "PG", "LPG", "PA", "LPA"];
        
        println!("\n");
        println!("Enter in the lipid species to calculate its accurate mass and adducts (e.g. PC-34:2):");
        println!("Options include: TG, CE, MG, DG, PC, LPC, PE, LPE, PS, LPS, PI, LPI, PG, LPG, PA, LPA");
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
            "DG" => dg_calc(&fa_carbon_no, &fa_desat_no),
            "MG" => mg_calc(&fa_carbon_no, &fa_desat_no),
            "CE" => ce_calc(&fa_carbon_no, &fa_desat_no),
            "PC" => pc_calc(&fa_carbon_no, &fa_desat_no),
            "PE" => pe_calc(&fa_carbon_no, &fa_desat_no),
            "PI" => pi_calc(&fa_carbon_no, &fa_desat_no),
            "PG" => pg_calc(&fa_carbon_no, &fa_desat_no),
            "PS" => ps_calc(&fa_carbon_no, &fa_desat_no),
            "PA" => pa_calc(&fa_carbon_no, &fa_desat_no),
            "LPC" => lpc_calc(&fa_carbon_no, &fa_desat_no),
            "LPE" => lpe_calc(&fa_carbon_no, &fa_desat_no),
            "LPI" => lpi_calc(&fa_carbon_no, &fa_desat_no),
            "LPG" => lpg_calc(&fa_carbon_no, &fa_desat_no),
            "LPS" => lps_calc(&fa_carbon_no, &fa_desat_no),
            "LPA" => lpa_calc(&fa_carbon_no, &fa_desat_no),
            _ => {
                println!("That lipid class hasn't been defined yet! Hold tight buckaroo!");
                continue;
            },
        
        };

        let positive_classes = vec!["TG", "CE", "MG", "DG"];
        let negative_classes = vec!["PC", "LPC", "PE", "LPE", "PS", "LPS", "PI", "LPI", "PG", "LPG", "PA", "LPA"];

        if positive_classes.contains(&lipid_class) {

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

        } else if negative_classes.contains(&lipid_class) {

            let monoiso_mass = &calc_masses[0];
            let h_minus_mass = &calc_masses[1];
            let ac_minus_mass = &calc_masses[2];
            let h_mass = &calc_masses[3];
            let na_mass = &calc_masses[4];
            let k_mass = &calc_masses[5];
            let nh4_mass = &calc_masses[6];

            println!("\n");
            println!("Monoisotopic mass: {monoiso_mass}");
            println!("[M-H]-: {h_minus_mass}");
            println!("[M+Ac-H]-: {ac_minus_mass}");
            println!("\n");
            println!("[M+H]+: {h_mass}");
            println!("[M+Na]+: {na_mass}");
            println!("[M+K]+: {k_mass}");
            println!("[M+NH4]+: {nh4_mass}");

            break;

        }


        
    };

        fn tg_calc(carb: &f64, desat: &f64) -> Vec<String> {
            let mut calculated_masses = Vec::new();

            let monoiso_mass: f64 = (CARBON * (carb + 3.0)) + (HYDROGEN * ((carb * 2.0) - (desat * 2.0) + 2.0)) + (OXYGEN * 6.0);

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

        fn dg_calc(carb: &f64, desat: &f64) -> Vec<String> {
            let mut calculated_masses = Vec::new();

            let monoiso_mass: f64 = (CARBON * (carb + 3.0)) + (HYDROGEN * ((carb * 2.0) - (desat * 2.0) + 4.0)) + (OXYGEN * 5.0);

            let dg_h_adduct = &monoiso_mass + HYDROGEN;
            let dg_na_adduct = &monoiso_mass + SODIUM;
            let dg_k_adduct = &monoiso_mass + POTASSIUM;
            let dg_nh4_adduct = &monoiso_mass + AMMONIUM;

            let monoiso_mass_form = format!("{monoiso_mass:.4}");
            let h_adduct_form = format!("{dg_h_adduct:.4}");
            let na_adduct_form = format!("{dg_na_adduct:.4}");
            let k_adduct_form = format!("{dg_k_adduct:.4}");
            let nh4_adduct_form = format!("{dg_nh4_adduct:.4}");

            calculated_masses.push(monoiso_mass_form);
            calculated_masses.push(h_adduct_form);
            calculated_masses.push(na_adduct_form);
            calculated_masses.push(k_adduct_form);
            calculated_masses.push(nh4_adduct_form);

            calculated_masses
        }

        fn mg_calc(carb: &f64, desat: &f64) -> Vec<String> {
            let mut calculated_masses = Vec::new();

            let monoiso_mass: f64 = (CARBON * (carb + 3.0)) + (HYDROGEN * ((carb * 2.0) - (desat * 2.0) + 6.0)) + (OXYGEN * 4.0);

            let mg_h_adduct = &monoiso_mass + HYDROGEN;
            let mg_na_adduct = &monoiso_mass + SODIUM;
            let mg_k_adduct = &monoiso_mass + POTASSIUM;
            let mg_nh4_adduct = &monoiso_mass + AMMONIUM;

            let monoiso_mass_form = format!("{monoiso_mass:.4}");
            let h_adduct_form = format!("{mg_h_adduct:.4}");
            let na_adduct_form = format!("{mg_na_adduct:.4}");
            let k_adduct_form = format!("{mg_k_adduct:.4}");
            let nh4_adduct_form = format!("{mg_nh4_adduct:.4}");

            calculated_masses.push(monoiso_mass_form);
            calculated_masses.push(h_adduct_form);
            calculated_masses.push(na_adduct_form);
            calculated_masses.push(k_adduct_form);
            calculated_masses.push(nh4_adduct_form);

            calculated_masses
        }

        fn ce_calc(carb: &f64, desat: &f64) -> Vec<String> {
            let mut calculated_masses = Vec::new();

            let monoiso_mass: f64 = ((CARBON * 27.0) + (HYDROGEN * 45.0) + (OXYGEN )) + (CARBON * carb) + (HYDROGEN * ((carb * 2.0) - (desat * 2.0) - 1.0)) + OXYGEN;

            let ce_h_adduct = &monoiso_mass + HYDROGEN;
            let ce_na_adduct = &monoiso_mass + SODIUM;
            let ce_k_adduct = &monoiso_mass + POTASSIUM;
            let ce_nh4_adduct = &monoiso_mass + AMMONIUM;

            let monoiso_mass_form = format!("{monoiso_mass:.4}");
            let h_adduct_form = format!("{ce_h_adduct:.4}");
            let na_adduct_form = format!("{ce_na_adduct:.4}");
            let k_adduct_form = format!("{ce_k_adduct:.4}");
            let nh4_adduct_form = format!("{ce_nh4_adduct:.4}");

            calculated_masses.push(monoiso_mass_form);
            calculated_masses.push(h_adduct_form);
            calculated_masses.push(na_adduct_form);
            calculated_masses.push(k_adduct_form);
            calculated_masses.push(nh4_adduct_form);

            calculated_masses
        }

        fn pc_calc(carb: &f64, desat: &f64) -> Vec<String> {
            let mut calculated_masses = Vec::new();

            let monoiso_mass: f64 = (CARBON * (carb + 8.0)) + (HYDROGEN * ((carb * 2.0) - (desat * 2.0) + 16.0)) + (OXYGEN * 8.0) + NITROGEN + PHOSPH;

            let pc_h_minus = &monoiso_mass - HYDROGEN;
            let pc_ac_minus = &monoiso_mass - HYDROGEN + ((2.0 * CARBON) + (2.0 * OXYGEN) + (3.0 * HYDROGEN));
            let pc_h_adduct = &monoiso_mass + HYDROGEN;
            let pc_na_adduct = &monoiso_mass + SODIUM;
            let pc_k_adduct = &monoiso_mass + POTASSIUM;
            let pc_nh4_adduct = &monoiso_mass + AMMONIUM;

            let monoiso_mass_form = format!("{monoiso_mass:.4}");
            let h_minus_form = format!("{pc_h_minus:.4}");
            let ac_minus_form = format!("{pc_ac_minus:.4}");
            let h_adduct_form = format!("{pc_h_adduct:.4}");
            let na_adduct_form = format!("{pc_na_adduct:.4}");
            let k_adduct_form = format!("{pc_k_adduct:.4}");
            let nh4_adduct_form = format!("{pc_nh4_adduct:.4}");

            calculated_masses.push(monoiso_mass_form);
            calculated_masses.push(h_minus_form);
            calculated_masses.push(ac_minus_form);
            calculated_masses.push(h_adduct_form);
            calculated_masses.push(na_adduct_form);
            calculated_masses.push(k_adduct_form);
            calculated_masses.push(nh4_adduct_form);

            calculated_masses
        }

        fn pe_calc(carb: &f64, desat: &f64) -> Vec<String> {
            let mut calculated_masses = Vec::new();

            let monoiso_mass: f64 = (CARBON * (carb + 5.0)) + (HYDROGEN * ((carb * 2.0) - (desat * 2.0) + 10.0)) + (OXYGEN * 8.0) + NITROGEN + PHOSPH;

            let pe_h_minus = &monoiso_mass - HYDROGEN;
            let pe_ac_minus = &monoiso_mass - HYDROGEN + ((2.0 * CARBON) + (2.0 * OXYGEN) + (3.0 * HYDROGEN));
            let pe_h_adduct = &monoiso_mass + HYDROGEN;
            let pe_na_adduct = &monoiso_mass + SODIUM;
            let pe_k_adduct = &monoiso_mass + POTASSIUM;
            let pe_nh4_adduct = &monoiso_mass + AMMONIUM;

            let monoiso_mass_form = format!("{monoiso_mass:.4}");
            let h_minus_form = format!("{pe_h_minus:.4}");
            let ac_minus_form = format!("{pe_ac_minus:.4}");
            let h_adduct_form = format!("{pe_h_adduct:.4}");
            let na_adduct_form = format!("{pe_na_adduct:.4}");
            let k_adduct_form = format!("{pe_k_adduct:.4}");
            let nh4_adduct_form = format!("{pe_nh4_adduct:.4}");

            calculated_masses.push(monoiso_mass_form);
            calculated_masses.push(h_minus_form);
            calculated_masses.push(ac_minus_form);
            calculated_masses.push(h_adduct_form);
            calculated_masses.push(na_adduct_form);
            calculated_masses.push(k_adduct_form);
            calculated_masses.push(nh4_adduct_form);

            calculated_masses
        }

        fn pi_calc(carb: &f64, desat: &f64) -> Vec<String> {
            let mut calculated_masses = Vec::new();

            let monoiso_mass: f64 = (CARBON * (carb + 9.0)) + (HYDROGEN * ((carb * 2.0) - (desat * 2.0) + 15.0)) + (OXYGEN * 13.0) + PHOSPH;

            let pi_h_minus = &monoiso_mass - HYDROGEN;
            let pi_ac_minus = &monoiso_mass - HYDROGEN + ((2.0 * CARBON) + (2.0 * OXYGEN) + (3.0 * HYDROGEN));
            let pi_h_adduct = &monoiso_mass + HYDROGEN;
            let pi_na_adduct = &monoiso_mass + SODIUM;
            let pi_k_adduct = &monoiso_mass + POTASSIUM;
            let pi_nh4_adduct = &monoiso_mass + AMMONIUM;

            let monoiso_mass_form = format!("{monoiso_mass:.4}");
            let h_minus_form = format!("{pi_h_minus:.4}");
            let ac_minus_form = format!("{pi_ac_minus:.4}");
            let h_adduct_form = format!("{pi_h_adduct:.4}");
            let na_adduct_form = format!("{pi_na_adduct:.4}");
            let k_adduct_form = format!("{pi_k_adduct:.4}");
            let nh4_adduct_form = format!("{pi_nh4_adduct:.4}");

            calculated_masses.push(monoiso_mass_form);
            calculated_masses.push(h_minus_form);
            calculated_masses.push(ac_minus_form);
            calculated_masses.push(h_adduct_form);
            calculated_masses.push(na_adduct_form);
            calculated_masses.push(k_adduct_form);
            calculated_masses.push(nh4_adduct_form);

            calculated_masses
        }

        fn pg_calc(carb: &f64, desat: &f64) -> Vec<String> {
            let mut calculated_masses = Vec::new();

            let monoiso_mass: f64 = (CARBON * (carb + 6.0)) + (HYDROGEN * ((carb * 2.0) - (desat * 2.0) + 11.0)) + (OXYGEN * 10.0) + PHOSPH;

            let pg_h_minus = &monoiso_mass - HYDROGEN;
            let pg_ac_minus = &monoiso_mass - HYDROGEN + ((2.0 * CARBON) + (2.0 * OXYGEN) + (3.0 * HYDROGEN));
            let pg_h_adduct = &monoiso_mass + HYDROGEN;
            let pg_na_adduct = &monoiso_mass + SODIUM;
            let pg_k_adduct = &monoiso_mass + POTASSIUM;
            let pg_nh4_adduct = &monoiso_mass + AMMONIUM;

            let monoiso_mass_form = format!("{monoiso_mass:.4}");
            let h_minus_form = format!("{pg_h_minus:.4}");
            let ac_minus_form = format!("{pg_ac_minus:.4}");
            let h_adduct_form = format!("{pg_h_adduct:.4}");
            let na_adduct_form = format!("{pg_na_adduct:.4}");
            let k_adduct_form = format!("{pg_k_adduct:.4}");
            let nh4_adduct_form = format!("{pg_nh4_adduct:.4}");

            calculated_masses.push(monoiso_mass_form);
            calculated_masses.push(h_minus_form);
            calculated_masses.push(ac_minus_form);
            calculated_masses.push(h_adduct_form);
            calculated_masses.push(na_adduct_form);
            calculated_masses.push(k_adduct_form);
            calculated_masses.push(nh4_adduct_form);

            calculated_masses
        }

        fn ps_calc(carb: &f64, desat: &f64) -> Vec<String> {
            let mut calculated_masses = Vec::new();

            let monoiso_mass: f64 = (CARBON * (carb + 6.0)) + (HYDROGEN * ((carb * 2.0) - (desat * 2.0) + 10.0)) + (OXYGEN * 10.0) + PHOSPH + NITROGEN;

            let ps_h_minus = &monoiso_mass - HYDROGEN;
            let ps_ac_minus = &monoiso_mass - HYDROGEN + ((2.0 * CARBON) + (2.0 * OXYGEN) + (3.0 * HYDROGEN));
            let ps_h_adduct = &monoiso_mass + HYDROGEN;
            let ps_na_adduct = &monoiso_mass + SODIUM;
            let ps_k_adduct = &monoiso_mass + POTASSIUM;
            let ps_nh4_adduct = &monoiso_mass + AMMONIUM;

            let monoiso_mass_form = format!("{monoiso_mass:.4}");
            let h_minus_form = format!("{ps_h_minus:.4}");
            let ac_minus_form = format!("{ps_ac_minus:.4}");
            let h_adduct_form = format!("{ps_h_adduct:.4}");
            let na_adduct_form = format!("{ps_na_adduct:.4}");
            let k_adduct_form = format!("{ps_k_adduct:.4}");
            let nh4_adduct_form = format!("{ps_nh4_adduct:.4}");

            calculated_masses.push(monoiso_mass_form);
            calculated_masses.push(h_minus_form);
            calculated_masses.push(ac_minus_form);
            calculated_masses.push(h_adduct_form);
            calculated_masses.push(na_adduct_form);
            calculated_masses.push(k_adduct_form);
            calculated_masses.push(nh4_adduct_form);

            calculated_masses
        }

        fn pa_calc(carb: &f64, desat: &f64) -> Vec<String> {
            let mut calculated_masses = Vec::new();

            let monoiso_mass: f64 = (CARBON * (carb + 3.0)) + (HYDROGEN * ((carb * 2.0) - (desat * 2.0) + 5.0)) + (OXYGEN * 8.0) + PHOSPH;

            let pa_h_minus = &monoiso_mass - HYDROGEN;
            let pa_ac_minus = &monoiso_mass - HYDROGEN + ((2.0 * CARBON) + (2.0 * OXYGEN) + (3.0 * HYDROGEN));
            let pa_h_adduct = &monoiso_mass + HYDROGEN;
            let pa_na_adduct = &monoiso_mass + SODIUM;
            let pa_k_adduct = &monoiso_mass + POTASSIUM;
            let pa_nh4_adduct = &monoiso_mass + AMMONIUM;

            let monoiso_mass_form = format!("{monoiso_mass:.4}");
            let h_minus_form = format!("{pa_h_minus:.4}");
            let ac_minus_form = format!("{pa_ac_minus:.4}");
            let h_adduct_form = format!("{pa_h_adduct:.4}");
            let na_adduct_form = format!("{pa_na_adduct:.4}");
            let k_adduct_form = format!("{pa_k_adduct:.4}");
            let nh4_adduct_form = format!("{pa_nh4_adduct:.4}");

            calculated_masses.push(monoiso_mass_form);
            calculated_masses.push(h_minus_form);
            calculated_masses.push(ac_minus_form);
            calculated_masses.push(h_adduct_form);
            calculated_masses.push(na_adduct_form);
            calculated_masses.push(k_adduct_form);
            calculated_masses.push(nh4_adduct_form);

            calculated_masses
        }

        fn lpc_calc(carb: &f64, desat: &f64) -> Vec<String> {
            let mut calculated_masses = Vec::new();

            let monoiso_mass: f64 = (CARBON * (carb + 8.0)) + (HYDROGEN * ((carb * 2.0) - (desat * 2.0) + 18.0)) + (OXYGEN * 7.0) + PHOSPH + NITROGEN;

            let lpc_h_minus = &monoiso_mass - HYDROGEN;
            let lpc_ac_minus = &monoiso_mass - HYDROGEN + ((2.0 * CARBON) + (2.0 * OXYGEN) + (3.0 * HYDROGEN));
            let lpc_h_adduct = &monoiso_mass + HYDROGEN;
            let lpc_na_adduct = &monoiso_mass + SODIUM;
            let lpc_k_adduct = &monoiso_mass + POTASSIUM;
            let lpc_nh4_adduct = &monoiso_mass + AMMONIUM;

            let monoiso_mass_form = format!("{monoiso_mass:.4}");
            let h_minus_form = format!("{lpc_h_minus:.4}");
            let ac_minus_form = format!("{lpc_ac_minus:.4}");
            let h_adduct_form = format!("{lpc_h_adduct:.4}");
            let na_adduct_form = format!("{lpc_na_adduct:.4}");
            let k_adduct_form = format!("{lpc_k_adduct:.4}");
            let nh4_adduct_form = format!("{lpc_nh4_adduct:.4}");

            calculated_masses.push(monoiso_mass_form);
            calculated_masses.push(h_minus_form);
            calculated_masses.push(ac_minus_form);
            calculated_masses.push(h_adduct_form);
            calculated_masses.push(na_adduct_form);
            calculated_masses.push(k_adduct_form);
            calculated_masses.push(nh4_adduct_form);

            calculated_masses
        }

        fn lpe_calc(carb: &f64, desat: &f64) -> Vec<String> {
            let mut calculated_masses = Vec::new();

            let monoiso_mass: f64 = (CARBON * (carb + 5.0)) + (HYDROGEN * ((carb * 2.0) - (desat * 2.0) + 12.0)) + (OXYGEN * 7.0) + PHOSPH + NITROGEN;

            let lpe_h_minus = &monoiso_mass - HYDROGEN;
            let lpe_ac_minus = &monoiso_mass - HYDROGEN + ((2.0 * CARBON) + (2.0 * OXYGEN) + (3.0 * HYDROGEN));
            let lpe_h_adduct = &monoiso_mass + HYDROGEN;
            let lpe_na_adduct = &monoiso_mass + SODIUM;
            let lpe_k_adduct = &monoiso_mass + POTASSIUM;
            let lpe_nh4_adduct = &monoiso_mass + AMMONIUM;

            let monoiso_mass_form = format!("{monoiso_mass:.4}");
            let h_minus_form = format!("{lpe_h_minus:.4}");
            let ac_minus_form = format!("{lpe_ac_minus:.4}");
            let h_adduct_form = format!("{lpe_h_adduct:.4}");
            let na_adduct_form = format!("{lpe_na_adduct:.4}");
            let k_adduct_form = format!("{lpe_k_adduct:.4}");
            let nh4_adduct_form = format!("{lpe_nh4_adduct:.4}");

            calculated_masses.push(monoiso_mass_form);
            calculated_masses.push(h_minus_form);
            calculated_masses.push(ac_minus_form);
            calculated_masses.push(h_adduct_form);
            calculated_masses.push(na_adduct_form);
            calculated_masses.push(k_adduct_form);
            calculated_masses.push(nh4_adduct_form);

            calculated_masses
        }

        fn lpi_calc(carb: &f64, desat: &f64) -> Vec<String> {
            let mut calculated_masses = Vec::new();

            let monoiso_mass: f64 = (CARBON * (carb + 9.0)) + (HYDROGEN * ((carb * 2.0) - (desat * 2.0) + 17.0)) + (OXYGEN * 12.0) + PHOSPH;

            let lpi_h_minus = &monoiso_mass - HYDROGEN;
            let lpi_ac_minus = &monoiso_mass - HYDROGEN + ((2.0 * CARBON) + (2.0 * OXYGEN) + (3.0 * HYDROGEN));
            let lpi_h_adduct = &monoiso_mass + HYDROGEN;
            let lpi_na_adduct = &monoiso_mass + SODIUM;
            let lpi_k_adduct = &monoiso_mass + POTASSIUM;
            let lpi_nh4_adduct = &monoiso_mass + AMMONIUM;

            let monoiso_mass_form = format!("{monoiso_mass:.4}");
            let h_minus_form = format!("{lpi_h_minus:.4}");
            let ac_minus_form = format!("{lpi_ac_minus:.4}");
            let h_adduct_form = format!("{lpi_h_adduct:.4}");
            let na_adduct_form = format!("{lpi_na_adduct:.4}");
            let k_adduct_form = format!("{lpi_k_adduct:.4}");
            let nh4_adduct_form = format!("{lpi_nh4_adduct:.4}");

            calculated_masses.push(monoiso_mass_form);
            calculated_masses.push(h_minus_form);
            calculated_masses.push(ac_minus_form);
            calculated_masses.push(h_adduct_form);
            calculated_masses.push(na_adduct_form);
            calculated_masses.push(k_adduct_form);
            calculated_masses.push(nh4_adduct_form);

            calculated_masses
        }

        fn lpg_calc(carb: &f64, desat: &f64) -> Vec<String> {
            let mut calculated_masses = Vec::new();

            let monoiso_mass: f64 = (CARBON * (carb + 6.0)) + (HYDROGEN * ((carb * 2.0) - (desat * 2.0) + 13.0)) + (OXYGEN * 9.0) + PHOSPH;

            let lpg_h_minus = &monoiso_mass - HYDROGEN;
            let lpg_ac_minus = &monoiso_mass - HYDROGEN + ((2.0 * CARBON) + (2.0 * OXYGEN) + (3.0 * HYDROGEN));
            let lpg_h_adduct = &monoiso_mass + HYDROGEN;
            let lpg_na_adduct = &monoiso_mass + SODIUM;
            let lpg_k_adduct = &monoiso_mass + POTASSIUM;
            let lpg_nh4_adduct = &monoiso_mass + AMMONIUM;

            let monoiso_mass_form = format!("{monoiso_mass:.4}");
            let h_minus_form = format!("{lpg_h_minus:.4}");
            let ac_minus_form = format!("{lpg_ac_minus:.4}");
            let h_adduct_form = format!("{lpg_h_adduct:.4}");
            let na_adduct_form = format!("{lpg_na_adduct:.4}");
            let k_adduct_form = format!("{lpg_k_adduct:.4}");
            let nh4_adduct_form = format!("{lpg_nh4_adduct:.4}");

            calculated_masses.push(monoiso_mass_form);
            calculated_masses.push(h_minus_form);
            calculated_masses.push(ac_minus_form);
            calculated_masses.push(h_adduct_form);
            calculated_masses.push(na_adduct_form);
            calculated_masses.push(k_adduct_form);
            calculated_masses.push(nh4_adduct_form);

            calculated_masses
        }

        fn lps_calc(carb: &f64, desat: &f64) -> Vec<String> {
            let mut calculated_masses = Vec::new();

            let monoiso_mass: f64 = (CARBON * (carb + 6.0)) + (HYDROGEN * ((carb * 2.0) - (desat * 2.0) + 12.0)) + (OXYGEN * 9.0) + PHOSPH + NITROGEN;

            let lps_h_minus = &monoiso_mass - HYDROGEN;
            let lps_ac_minus = &monoiso_mass - HYDROGEN + ((2.0 * CARBON) + (2.0 * OXYGEN) + (3.0 * HYDROGEN));
            let lps_h_adduct = &monoiso_mass + HYDROGEN;
            let lps_na_adduct = &monoiso_mass + SODIUM;
            let lps_k_adduct = &monoiso_mass + POTASSIUM;
            let lps_nh4_adduct = &monoiso_mass + AMMONIUM;

            let monoiso_mass_form = format!("{monoiso_mass:.4}");
            let h_minus_form = format!("{lps_h_minus:.4}");
            let ac_minus_form = format!("{lps_ac_minus:.4}");
            let h_adduct_form = format!("{lps_h_adduct:.4}");
            let na_adduct_form = format!("{lps_na_adduct:.4}");
            let k_adduct_form = format!("{lps_k_adduct:.4}");
            let nh4_adduct_form = format!("{lps_nh4_adduct:.4}");

            calculated_masses.push(monoiso_mass_form);
            calculated_masses.push(h_minus_form);
            calculated_masses.push(ac_minus_form);
            calculated_masses.push(h_adduct_form);
            calculated_masses.push(na_adduct_form);
            calculated_masses.push(k_adduct_form);
            calculated_masses.push(nh4_adduct_form);

            calculated_masses
        }

        fn lpa_calc(carb: &f64, desat: &f64) -> Vec<String> {
            let mut calculated_masses = Vec::new();

            let monoiso_mass: f64 = (CARBON * (carb + 3.0)) + (HYDROGEN * ((carb * 2.0) - (desat * 2.0) + 7.0)) + (OXYGEN * 7.0) + PHOSPH;

            let lpa_h_minus = &monoiso_mass - HYDROGEN;
            let lpa_ac_minus = &monoiso_mass - HYDROGEN + ((2.0 * CARBON) + (2.0 * OXYGEN) + (3.0 * HYDROGEN));
            let lpa_h_adduct = &monoiso_mass + HYDROGEN;
            let lpa_na_adduct = &monoiso_mass + SODIUM;
            let lpa_k_adduct = &monoiso_mass + POTASSIUM;
            let lpa_nh4_adduct = &monoiso_mass + AMMONIUM;

            let monoiso_mass_form = format!("{monoiso_mass:.4}");
            let h_minus_form = format!("{lpa_h_minus:.4}");
            let ac_minus_form = format!("{lpa_ac_minus:.4}");
            let h_adduct_form = format!("{lpa_h_adduct:.4}");
            let na_adduct_form = format!("{lpa_na_adduct:.4}");
            let k_adduct_form = format!("{lpa_k_adduct:.4}");
            let nh4_adduct_form = format!("{lpa_nh4_adduct:.4}");

            calculated_masses.push(monoiso_mass_form);
            calculated_masses.push(h_minus_form);
            calculated_masses.push(ac_minus_form);
            calculated_masses.push(h_adduct_form);
            calculated_masses.push(na_adduct_form);
            calculated_masses.push(k_adduct_form);
            calculated_masses.push(nh4_adduct_form);

            calculated_masses
        }
}