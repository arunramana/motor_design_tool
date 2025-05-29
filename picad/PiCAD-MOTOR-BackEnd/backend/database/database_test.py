from pymongo import MongoClient
import json
from bson import ObjectId
import gridfs

import warnings
warnings.filterwarnings("ignore")

# Connect to MongoDB
#connection_str = 'mongodb://Admin28:sK_25171@moto.cluster-crwa8q8k0lvn.ap-south-1.docdb.amazonaws.com:27017/?tls=true&tlsCAFile=global-bundle.pem&replicaSet=rs0&readPreference=secondaryPreferred&retryWrites=false'

connection_str = "mongodb://65.0.168.91:27017/"

client = MongoClient(connection_str)

db_users = client['users']
db_vehicles = client['vehicles']
db_applications = client['applications']
db_motors = client['motors']

db_emails = client['emails']




db_downloads = client['downloads']



# Create collections (tables)
users_collection = db_users['users']
vehicles_collection = db_vehicles['vehicles']
applications_collection = db_applications['applications']
motors_collection = db_motors['motors']

emails_collection = db_emails['emails']

downloads_collection = db_downloads['downloads']


fs_downloads = gridfs.GridFS(db_downloads)

# Define schema for collections
users_schema = {
    "user_name": str,
    "email": str,
    "password": str
}

vehicles_schema = {
    "user_id": ObjectId,
    "vehicle_name": str,
    "input": object,
    "output": object,
    "error": str,
    "desc": str
} 

#for non ev application
application_schema = {
    "user_id": ObjectId,
    "application_name": str,
    "input": object,
    "output": object,
    "error": str,
    "desc": str
} 

motors_schema = {
    "user_id": ObjectId,
    "vehicle_name": str,
    "motor_name": str,
    "input": object,
    "output": object,
    "error": str,
    "desc": str,
    "status":str,
    "is_ev":bool
}


emails_schema={
    "email": str
}


student_emails = ['mc010524001@code.iitm.ac.in',
 'mc010524002@code.iitm.ac.in',
 'mc010524003@code.iitm.ac.in',
 'mc010524004@code.iitm.ac.in',
 'mc010524005@code.iitm.ac.in',
 'mc010524006@code.iitm.ac.in',
 'mc010524007@code.iitm.ac.in',
 'mc010524008@code.iitm.ac.in',
 'mc010524009@code.iitm.ac.in',
 'mc010524010@code.iitm.ac.in',
 'mc010524011@code.iitm.ac.in',
 'mc010524012@code.iitm.ac.in',
 'mc010524013@code.iitm.ac.in',
 'mc010524014@code.iitm.ac.in',
 'mc010524015@code.iitm.ac.in',
 'mc010524016@code.iitm.ac.in',
 'mc010524017@code.iitm.ac.in',
 'mc010524018@code.iitm.ac.in',
 'mc010524019@code.iitm.ac.in',
 'mc010524020@code.iitm.ac.in',
 'mc010524021@code.iitm.ac.in',
 'mc010524022@code.iitm.ac.in',
 'mc010524023@code.iitm.ac.in',
 'mc010524024@code.iitm.ac.in',
 'mc010524025@code.iitm.ac.in',
 'mc010524026@code.iitm.ac.in',
 'mc010524027@code.iitm.ac.in',
 'mc010524028@code.iitm.ac.in',
 'mc010524029@code.iitm.ac.in',
 'mc010524030@code.iitm.ac.in',
 'mc010524031@code.iitm.ac.in',
 'mc010524032@code.iitm.ac.in',
 'mc010524033@code.iitm.ac.in',
 'mc010524034@code.iitm.ac.in',
 'mc010524035@code.iitm.ac.in',
 'mc010524036@code.iitm.ac.in',
 'mc010524037@code.iitm.ac.in',
 'mc010524038@code.iitm.ac.in',
 'mc010524039@code.iitm.ac.in',
 'mc010524040@code.iitm.ac.in',
 'mc010524041@code.iitm.ac.in',
 'mc010524042@code.iitm.ac.in',
 'mc010524043@code.iitm.ac.in',
 'testmoco@email.com']

emails_pilabz = [
    "jitin.g@pilabz.in",
    "beeulah.m@pilabz.in",
    "inbaraj.k@pilabz.in",
    "antony.j@pilabz.in",
    "aravindhan.a@pilabz.in",
    "balamuralikrishna.d@pilabz.in",
    "senthilkumar.k@pilabz.in",
    "yuvaraja.t@pilabz.in",
    "ponclinton.a@pilabz.in",
    "oscaranandh.s@pilabz.in",
    "niraipandi.p@pilabz.in",
    "vigneshwar.vw@pilabz.in",
    "supreeth.ss@pilabz.in",
    "venkatesh.sankar@pilabz.in",
    "sharikh.m@pilabz.in",
    "janani.ra@pilabz.in",
    "jagannath.vm@pilabz.in",
    "albin.saji@pilabz.in",
    "manirathnam.ts@pilabz.in",
    "aswanth.mohanan@pilabz.in",
    "pradeep.ak@pilabz.in",
    "balaji.tkk@pilabz.in",
    "mustafa.alam@pilabz.in",
    "chinmaya.vijay@pilabz.in",
    "tangirala.s@pilabz.in",
    "kannan@pilabz.in",
    "hemamalini.j@pilabz.in",
    "guruvignesh@pilabz.in"
]


#display emails
for x in emails_collection.find({}):
    print(x)


def create_email(email):
    try:
        # Insert email into emails collection
        email_id = emails_collection.insert_one({"email": email}).inserted_id
        return email_id, 200
    except Exception as e:
        return  str(e), 500
    

#un comment to insert data (run only once!!!!! at the beginning!!!!!)

'''
#insert student emails
for e in student_emails:
    create_email(e)
    
    
#insert pilabz email
for e in emails_pilabz:
    create_email(e)
'''    
    


def delete_email(email):
    try:
        # Delete email from emails collection
        result = emails_collection.delete_one({"email": email})
        return result.deleted_count, 200
    
    except Exception as e:
        return str(e), 500

def get_email(email):
    try:
        # Find email in emails collection
        email_data = emails_collection.find_one({"email": email})
        if email_data:
            return email_data["email"], 200
        else:
            return "Email not found", 404
    except Exception as e:
        return str(e), 500

def update_email(current_email, new_email):
    try:
        # Update email in emails collection
        result = emails_collection.update_one({"email": current_email}, {"$set": {"email": new_email}})
        return result.modified_count, 200
    except Exception as e:
        return str(e), 500

def create_user(username, email, password):
    try:
        user_name = username
        email = email
        password = password
        
        # Insert user into users collection
        if(users_collection.find_one({"email": email, "password": password})):
            return "User Already Exists", 500
        
        user_id = users_collection.insert_one({"user_name": user_name, "email": email, "password": password}).inserted_id
        return [str(user_id), str(email), str(user_name)], 200
    except Exception as e:
        return str(e), 500

def delete_user(user_id):
    try:
        # Delete user
        result_user = users_collection.delete_one({"_id": ObjectId(user_id)})
    
        # Get the vehicle names to delete related motors
        vehicle_names = [vehicle['vehicle_name'] for vehicle in vehicles_collection.find({"user_id": ObjectId(user_id)}, {"vehicle_name": 1})]
    
        # Delete motors related to the vehicles
        result_motors = motors_collection.delete_many({"vehicle_name": {"$in": vehicle_names}})
    
        # Delete vehicles related to the user
        result_vehicles = vehicles_collection.delete_many({"user_id": ObjectId(user_id)})
    
        return json.dumps({"deleted_user_count": result_user.deleted_count,
                           "deleted_vehicle_count": result_vehicles.deleted_count,
                           "deleted_motor_count": result_motors.deleted_count})
    except Exception as e:
        return str(e), 500

def create_vehicle(user_id, vehicle_name, input_data, output_data, error, desc, is_demo=False):
    try:
        # Insert vehicle into vehicles collection
        vehicle_id = vehicles_collection.insert_one({"user_id": ObjectId(user_id), "vehicle_name": vehicle_name, "input": input_data, "output": output_data, "error": error, "desc": desc}).inserted_id
    
        return str(vehicle_id), 200
    
    except Exception as e:
        return str(e), 500

def update_vehicle(user_id, vehicle_name, input_data, output_data, error_message):
    try:
        # Update vehicle in vehicles collection
        vehicles_collection.update_one({"user_id": ObjectId(user_id), "vehicle_name": vehicle_name}, {"$set": {"input": input_data, "output": output_data, "error": error_message}})
        
        # Delete all motors related to this vehicle
        motors_collection.delete_many({"user_id": ObjectId(user_id), "vehicle_name": vehicle_name})
        
        return "Vehicle updated successfully", 200
    except Exception as e:
        return str(e), 500


def delete_vehicle(user_id, vehicle_name):
    try:
        # Delete the vehicle
        result = vehicles_collection.delete_one({"vehicle_name": vehicle_name, "user_id": ObjectId(user_id)})
        
        # Check if the vehicle was found and deleted
        if result.deleted_count == 1:
            # Delete all motors related to the deleted vehicle
            motors_collection.delete_many({"user_id": ObjectId(user_id), "vehicle_name": vehicle_name})
            
            return "Vehicle deleted successfully", 200
        else:
            return "Vehicle not found", 404
    except Exception as e:
        return str(e), 500
        
        
def get_vehicle_data(user_id, vehicle_name):
    try:
        # Find the vehicle details for the given user_id and vehicle_name
        vehicle = vehicles_collection.find_one({"user_id": ObjectId(user_id), "vehicle_name": vehicle_name})
        
        # Check if vehicle exists
        if vehicle:
            # Extract vehicle details
            vehicle_data = {
                "vehicle_name": vehicle["vehicle_name"],
                "input": vehicle["input"],
                "output": vehicle["output"],
                "error": vehicle["error"],
                "desc": vehicle["desc"]
            }
            return vehicle_data, 200
        else:
            return "Vehicle not found", 404
    except Exception as e:
        return str(e), 500

        
def get_all_vehicles(user_id):
    try:
        # Find all vehicles for the given user_id
        vehicles = vehicles_collection.find({"user_id": ObjectId(user_id)}, {"vehicle_name": 1, "desc": 1})
        
        # Extract vehicle names and descriptions
        vehicle_data = [{"vehicle_name": vehicle["vehicle_name"], "desc": vehicle["desc"]} for vehicle in vehicles]
        
        return vehicle_data, 200
    except Exception as e:
        return str(e), 500
        
        
def is_existing_vehicle(user_id, vehicle_name):
    vehicle = vehicles_collection.find_one({"user_id": ObjectId(user_id), "vehicle_name": vehicle_name})
    if vehicle:
        return True
    else:
        return False

def create_vehicle_version(user_id, vehicle_name):
    try:
        # Check if the vehicle already exists
        existing_vehicle = vehicles_collection.find_one({"user_id": ObjectId(user_id), "vehicle_name": vehicle_name})
        
        if existing_vehicle:
            i = 2
            desc = existing_vehicle["desc"]
            
            while existing_vehicle:
                # Generate version name
                version_name = f"{vehicle_name}_version_{i}"
                
                # Check if the version name already exists
                existing_vehicle = vehicles_collection.find_one({"user_id": ObjectId(user_id), "vehicle_name": version_name})
                
                if not existing_vehicle:
                    # Update vehicle_name in the database
                    #vehicles_collection.update_one({"user_id": ObjectId(user_id), "vehicle_name": vehicle_name}, {"$set": {"vehicle_name": version_name}})
                    #return version_name
                    break
                else:
                    desc = existing_vehicle["desc"]
                    i += 1
            
            # create vehicle_name in the database
            create_vehicle(user_id, f"{vehicle_name}_version_{i}", [], [], "none", desc)
            return f"{vehicle_name}_version_{i}", 200
        else:
            return "Vehicle not found", 404
    except Exception as e:
        return str(e), 500


def create_application_version(user_id, application_name):
    try:
        # Check if the application already exists
        existing_application = applications_collection.find_one({"user_id": ObjectId(user_id), "application_name": application_name})
        
        if existing_application:
            i = 2
            desc = existing_application["desc"]
            
            while existing_application:
                # Generate version name
                version_name = f"{application_name}_version_{i}"
                
                # Check if the version name already exists
                existing_application = applications_collection.find_one({"user_id": ObjectId(user_id), "application_name": version_name})
                
                if not existing_application:
                    # Update application_name in the database
                    #applications_collection.update_one({"user_id": ObjectId(user_id), "application_name": application_name}, {"$set": {"application_name": version_name}})
                    #return version_name
                    break
                else:
                    desc = existing_application["desc"]
                    i += 1
            
            # create application_name in the database
            create_application(user_id, f"{application_name}_version_{i}", [], [], "none", desc)
            return f"{application_name}_version_{i}", 200
        else:
            return "application not found", 404
    except Exception as e:
        return str(e), 500

def create_motor(user_id, vehicle_name, motor_name, input_data, output_data, error, desc, status, is_ev=True):
    try:
        # Insert motor into motors collection
        motor_id = motors_collection.insert_one({"user_id": ObjectId(user_id), "vehicle_name": vehicle_name, "motor_name": motor_name, "input": input_data, "output": output_data, "error": error, "desc": desc, "status": status, "is_ev": is_ev}).inserted_id
    
        return str(motor_id), 200
    
    except Exception as e:
        return str(e), 500
        
        
def update_motor(user_id, vehicle_name, motor_name, input_data, output_data, error,status):
    try:
        # Update motor in motors collection
        result = motors_collection.update_one(
            {"user_id": ObjectId(user_id), "vehicle_name": vehicle_name, "motor_name": motor_name},
            {"$set": {"input": input_data, "output": output_data, "error": error, "status": status}}
        )
        
        if result.matched_count == 1:
            return "Motor updated successfully", 200
        else:
            return "Motor not found", 404
        
    except Exception as e:
        return str(e), 500


def delete_motor(user_id,vehicle_name,motor_name):
    try:
        result = motors_collection.delete_one({"vehicle_name": vehicle_name, "user_id": ObjectId(user_id), "motor_name": motor_name})
        return result.deleted_count, 200
    
    except Exception as e:
        return str(e), 500
        

def get_motor_data(user_id, vehicle_name, motor_name):
    try:
        # Find the motor details for the given user_id, vehicle_name, and motor_name
        motor = motors_collection.find_one({"user_id": ObjectId(user_id), "vehicle_name": vehicle_name, "motor_name": motor_name})
        
        # Check if motor exists
        if motor:
            # Extract motor details
            motor_data = {
                "motor_name": motor["motor_name"],
                "input": motor["input"],
                "output": motor["output"],
                "error": motor["error"],
                "desc": motor["desc"]
            }
            return motor_data, 200
        else:
            return "Motor not found", 404
    except Exception as e:
        return str(e), 500


def get_all_motors(user_id, vehicle_name):
    try:
        # Find all motors for the given user_id and vehicle_name
        motors = motors_collection.find({"user_id": ObjectId(user_id), "vehicle_name": vehicle_name}, {"motor_name": 1, "desc": 1})
        
        # Extract motor names and descriptions
        motor_data = [{"motor_name": motor["motor_name"], "desc": motor["desc"]} for motor in motors]
        
        return motor_data, 200
    except Exception as e:
        return str(e), 500
        
def delete_all_motors(user_id, vehicle_name):
    try:
        # Delete all motors related to the given user_id and vehicle_name
        result = motors_collection.delete_many({"user_id": ObjectId(user_id), "vehicle_name": vehicle_name})
        
        return result.deleted_count, 200
    except Exception as e:
        return str(e), 500
        
def create_motor_version(user_id, vehicle_name, motor_name):
    try:
        # Check if the motor already exists
        existing_motor = motors_collection.find_one({"user_id": ObjectId(user_id), "vehicle_name": vehicle_name, "motor_name": motor_name})
        
        if existing_motor:
            i = 2
            desc = existing_motor["desc"]
            inp = existing_motor["input"]
            op = existing_motor["output"]
            error = existing_motor["error"]
            status = existing_motor["status"]
            is_ev = existing_motor["is_ev"]
            
            while existing_motor:
                # Generate version name
                version_name = f"{motor_name}_version_{i}"
                
                # Check if the version name already exists
                existing_motor = motors_collection.find_one({"user_id": ObjectId(user_id), "vehicle_name": vehicle_name, "motor_name": version_name})
                
                if not existing_motor:
                    # Update motor_name in the database
                    #motors_collection.update_one({"user_id": ObjectId(user_id), "motor_name": motor_name}, {"$set": {"motor_name": version_name}})
                    #return version_name
                    break
                else:
                    desc = existing_motor["desc"]
                    inp = existing_motor["input"]
                    op = existing_motor["output"]
                    error = existing_motor["error"]
                    status = existing_motor["status"]
                    is_ev = existing_motor["is_ev"]
                    i += 1
            
            # create motor_name in the database
            
            create_motor(user_id, vehicle_name, f"{motor_name}_version_{i}", inp, op, error, desc, status, is_ev)
            
            
            return f"{motor_name}_version_{i}", 200
        else:
            return "motor not found", 404
    except Exception as e:
        return str(e), 500 

        
def is_existing_motor(user_id, vehicle_name, motor_name):
    motor = motors_collection.find_one({"user_id": ObjectId(user_id), "vehicle_name": vehicle_name, "motor_name": motor_name})
    if motor:
        return True
    else:
        return False



def login(email, password):
    try:
        # Search for the user in the database
        user = users_collection.find_one({"email": email, "password": password})
        
        if user:
            # User found, return the user ID
            return [str(user['_id']), str(user['email']), str(user['user_name'])], 200
        else:
            # User not found or invalid credentials
            return "Invalid email or password", 401
    except Exception as e:
        return str(e), 500
    
def store_zip_file(user_id, vehicle_name, motor_name, file_path):
    try:
        # Open the zip file and store it in GridFS within the downloads database
        with open(file_path, 'rb') as file_data:
            file_id = fs_downloads.put(file_data, filename=file_path)

        # Insert file metadata into the download_collection
        download_id = downloads_collection.insert_one({
            "user_id": ObjectId(user_id),
            "vehicle_name": vehicle_name,
            "motor_name": motor_name,
            "file_id": file_id,
            "filename": file_path
        }).inserted_id

        return str(download_id), 200
    except Exception as e:
        return str(e), 500  
        

# Function to retrieve zip file from downloads collection
def retrieve_zip_file_from_downloads(user_id, vehicle_name, motor_name):
    try:
        # Find the file metadata
        file_meta = downloads_collection.find_one({
            "user_id": ObjectId(user_id),
            "vehicle_name": vehicle_name,
            "motor_name": motor_name
        })

        if not file_meta:
            return json.dumps({"error": "File metadata not found"}), 404

        file_id = file_meta['file_id']
        filename = file_meta['filename']

        # Retrieve the file from GridFS
        file_data = fs_downloads.get(file_id)
        file_content = file_data.read()

       

        return file_content, 200
    except Exception as e:
        return str(e), 500


def create_application(user_id, application_name, input_data, output_data, error, desc):
    try:
        # Insert application into applications collection
        application_id = applications_collection.insert_one({
            "user_id": ObjectId(user_id),
            "application_name": application_name,
            "input": input_data,
            "output": output_data,
            "error": error,
            "desc": desc,
        }).inserted_id
    
        return str(application_id), 200
    except Exception as e:
        return str(e), 500

def update_application(user_id, application_name, input_data, output_data, error_message):
    try:
        # Update application in applications collection
        result = applications_collection.update_one(
            {"user_id": ObjectId(user_id), "application_name": application_name},
            {"$set": {"input": input_data, "output": output_data, "error": error_message}}
        )
        
        if result.matched_count == 1:
            return "Application updated successfully", 200
        else:
            return "Application not found", 404
    except Exception as e:
        return str(e), 500

def delete_application(user_id, application_name):
    try:
        # Delete the application
        result = applications_collection.delete_one({
            "user_id": ObjectId(user_id),
            "application_name": application_name
        })
        
        if result.deleted_count == 1:
            return "Application deleted successfully", 200
        else:
            return "Application not found", 404
    except Exception as e:
        return str(e), 500

def get_application_data(user_id, application_name):
    try:
        # Find the application details for the given user_id and application_name
        application = applications_collection.find_one({
            "user_id": ObjectId(user_id),
            "application_name": application_name
        })
        
        if application:
            # Extract application details
            application_data = {
                "application_name": application["application_name"],
                "input": application["input"],
                "output": application["output"],
                "error": application["error"],
                "desc": application["desc"]
            }
            return application_data, 200
        else:
            return "Application not found", 404
    except Exception as e:
        return str(e), 500

def get_all_applications(user_id):
    try:
        # Find all applications for the given user_id
        applications = applications_collection.find(
            {"user_id": ObjectId(user_id)},
            {"application_name": 1, "desc": 1}
        )
        
        # Extract application names and descriptions
        application_data = [{"application_name": application["application_name"], "desc": application["desc"]} for application in applications]
        
        return application_data, 200
    except Exception as e:
        return str(e), 500

def is_existing_application(user_id, application_name):
    try:
        application = applications_collection.find_one({
            "user_id": ObjectId(user_id),
            "application_name": application_name
        })
        return application is not None
    except Exception as e:
        return False
