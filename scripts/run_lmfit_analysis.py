import sys
import logging
from pathlib import Path
import lmfit

# Ensure project root is in sys.path
sys.path.append(str(Path(__file__).resolve().parent.parent))

from unified_optimizer import config
from unified_optimizer.data_manager import DataManager
from unified_optimizer.optimizer_lmfit import run_lmfit_2d

def main():
    logging.basicConfig(level=logging.INFO, format='%(asctime)s [%(levelname)s] %(message)s')
    
    data_dir = Path(config.DATA_DIR)
    manager = DataManager(data_dir)
    
    # We will run on the best dataset '356att'
    dataset = "356att"
    logging.info(f"Loading data for {dataset}...")
    data_dict = manager.get_data_for_dataset(dataset)
    
    if not data_dict:
        logging.error(f"Failed to load dataset {dataset}.")
        return

    result, mini = run_lmfit_2d(data_dict, dataset_name=dataset)
    
    report_text = lmfit.fit_report(result)
    print("\n" + "="*50)
    print("LMFIT REPORT:")
    print("="*50)
    print(report_text)
    
    try:
        logging.info("Calculating strict confidence intervals (this may take a moment)...")
        ci = lmfit.conf_interval(mini, result)
        print("\n" + "="*50)
        print("CONFIDENCE INTERVALS:")
        print("="*50)
        lmfit.printfuncs.report_ci(ci)
        
        ci_report = lmfit.printfuncs.ci_report(ci)
        report_text += "\n\n" + "="*50 + "\nCONFIDENCE INTERVALS\n" + "="*50 + "\n" + ci_report
    except Exception as e:
        logging.warning(f"Failed to compute confidence intervals: {e}")
        
    results_dir = Path(config.BASE_DIR) / "results"
    results_dir.mkdir(exist_ok=True)
    report_file = results_dir / "lmfit_report.txt"
    
    with open(report_file, "w", encoding="utf-8") as f:
        f.write(report_text)
        
    logging.info(f"Report saved to {report_file}")

if __name__ == "__main__":
    main()
