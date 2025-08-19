use std::io;
mod lipid_calculator;

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


