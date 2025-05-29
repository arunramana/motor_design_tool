from pymongo import MongoClient, ASCENDING
from datetime import datetime
import random
from datetime import datetime, timedelta
import threading
import time
# MongoDB connection settings for local testing


# Connect to MongoDB
try:
    #client = MongoClient('localhost', 27017)
    #connection_str = 'mongodb://Admin28:sK_25171@moto.cluster-crwa8q8k0lvn.ap-south-1.docdb.amazonaws.com:27017/?tls=true&tlsCAFile=backend/database/global-bundle.pem&replicaSet=rs0&readPreference=secondaryPreferred&retryWrites=false'
    
    connection_str = "mongodb://65.0.168.91:27017/"
    # connection_str = "mongodb://localhost:27017/"
    
    client = MongoClient(connection_str)

    client.admin.command('ping')  # Verify connection to the server
    #print("MongoDB connection established successfully.")
except Exception as e:
    print("Failed to connect to MongoDB:", e)
    raise e

# Define the User database and collection
db_users = client['users']
users_collection = db_users['users']



# Define the OTP database and collection
db_otp = client['otp']
otp_collection = db_otp['otp']

# Create a TTL index on the 'timestamp' field
# This will automatically delete documents 10 minutes (600 seconds) after the 'timestamp'
otp_collection.create_index([("timestamp", ASCENDING)], expireAfterSeconds=600)

# Function to generate a 6-digit OTP
def generate_otp():
    return random.randint(100000, 999999)

# Function to check if email exists in the users collection
def email_exists(email):
    return users_collection.find_one({"email": email}) is not None


def email_exists_in_otp(email):
    try:
        return otp_collection.find_one({"email": email}) is not None
    except Exception as e:
        return False





def store_otp1(email, otp):
    try:
        if otp_collection.find_one({"email": email}):
            timestamp = datetime.utcnow()
            result = otp_collection.update_one(
                {"email": email},
                {"$set": {"otp": otp, "timestamp": timestamp}}
            )
        else:
            otp_data = {
                "email": email,
                "otp": otp,
                "timestamp": datetime.utcnow()
            }
            result = otp_collection.insert_one(otp_data)
        return "OTP stored successfully", 200
    except Exception as e:
        return str(e), 500


# Function to store the OTP along with email and timestamp
def store_otp(email, otp):
    try:
        if not email_exists(email):
            return "Email not found in user database", 404
        
        
        if(otp_collection.find_one({"email": email}) is not None):
            
            #update otp for resend otp
            
            otp = otp
            timestamp = datetime.utcnow()
            otp_data = {
                "email": email,
                "otp": otp,
                "timestamp": timestamp
            }
            
            result = otp_collection.update_one({"email": email}, {"$set": {"otp": otp, "timestamp": timestamp}})
        
        else:
        
            #create otp
            
            otp = otp
            timestamp = datetime.utcnow()
            otp_data = {
                "email": email,
                "otp": otp,
                "timestamp": timestamp
            }
            result = otp_collection.insert_one(otp_data)
            
        return "OTP stored successfully", 200
    except Exception as e:
        return str(e), 500
    
def verify_otp1(email, otp):
    try:
        otp_data = otp_collection.find_one({"email": email, "otp": otp})
        if not otp_data:
            return "Invalid OTP or email", 400
        
        current_time = datetime.utcnow()
        otp_timestamp = otp_data["timestamp"]
        elapsed_time = (current_time - otp_timestamp).total_seconds()

        if elapsed_time > 600:  # OTP is valid for 10 minutes
            return "OTP expired", 400
        
        return  "OTP verified successfully", 200
    except Exception as e:
        return str(e), 500

# Function to verify the OTP
def verify_otp(email, otp):
    try:
        otp_data = otp_collection.find_one({"email": email, "otp": otp})
        if not otp_data:
            return "Invalid OTP or email", 400
        
        current_time = datetime.utcnow()
        otp_timestamp = otp_data["timestamp"]
        elapsed_time = (current_time - otp_timestamp).total_seconds()

        if elapsed_time > 600:  # OTP is valid for 10 minutes
            return "OTP expired", 400
        
        return  "OTP verified successfully", 200
    except Exception as e:
        return str(e), 500

# Function to update the password
def update_password(email, new_password):
    try:
        result = users_collection.update_one({"email": email}, {"$set": {"password": new_password}})
        if result.matched_count == 1:
            return "Password updated successfully", 200
        else:
            return "Failed to update password", 500
    except Exception as e:
        return  str(e), 500

# Example usage
def main():
    try:
        data=users_collection.find()
        for datam in data:
            print(datam)
        # Generate and store OTP
       # email = "bala@gmail.com"
        #store_response = store_otp(email)
        #print(store_response)

        #if "error" in store_response:
         #   return

        # Suppose the user received the OTP and enters it for verification
       # user_otp = store_response['otp']  # In practice, this would be input by the user

        # Verify OTP
       # verify_response = verify_otp(email, user_otp)
       # print(verify_response)

      #  if "error" in verify_response:
      #      return
        
        # Update password after OTP verification
     #   new_password = "newpassword123"  # In practice, this would be input by the user
      #  update_response = update_password(email, new_password)
       # print(update_response)
        

       
        
    #    client.close()
    except Exception as e:
        print("An error occurred:", e)
    
'''
if __name__ == "__main__":
    main()
    
'''